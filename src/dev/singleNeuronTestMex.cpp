#include "mex.h"
#include "matFunc.h"
#include "time.h"
#include <cstring>
#include "channels.h"
#include "NeuronLibrary.h"
#include "NeuronStruct.h"
#include "singleRK4.h"
#include "singleRK2.h"
#include "singleBilinear.h"
#include "typedefs.h"

using std::cout;
using std::endl;
using std::vector;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // prhs[] ~ lib_file, para_file, ith, rE, rI, model, run_t, ignore_t, vinit
    char *lib_file,*para_file;
    mxArray *para, *tmp;
    double cpu_t_sim, cpu_t_bilinear, cpu_t_linear;
    clockid_t clk_id = CLOCK_PROCESS_CPUTIME_ID;
    struct timespec tpS, tpE;
    unsigned int ith,i,j,k,model,nt,rk;
    double rE, rI, run_t, ignore_t, tau_ed, tau_er, tau_id, tau_ir;
    double gNa, vNa, gK, vK, gLeak, vLeak, vT, vinit, vE, vI, vRest, DeltaT, S;
    double tsp_sim, tsp_bi, tsp_li;
    MATFile *matFile;
    NeuroLib neuroLib;
    Neuron neuron;
    vector<double> fE,fI,fStrength;
    vector<double> simV, gE, gI, m, n, h, biV, liV;
    bool cutoff;
    unsigned int seed;
	//unsigned int seed = static_cast<unsigned int>(std::time(NULL));

    size_t lib_file_l = mxGetN(prhs[0])+1;
    lib_file = (char *) mxMalloc(lib_file_l);
    mxGetString(prhs[0],lib_file,(mwSize) lib_file_l);
    neuroLib.readLib(lib_file);

    size_t para_file_l = mxGetN(prhs[1])+1;
    para_file = (char *) mxMalloc(para_file_l);
    mxGetString(prhs[1],para_file,(mwSize) para_file_l);

    ith = static_cast<unsigned int>(mxGetScalar(prhs[2]))-1;
    rE = mxGetScalar(prhs[3]);
    rI = mxGetScalar(prhs[4]);
    model = static_cast<unsigned int>(mxGetScalar(prhs[5]));
    run_t = mxGetScalar(prhs[6]);
    ignore_t = mxGetScalar(prhs[7]);
    vinit = mxGetScalar(prhs[8]);
    rk = mxGetScalar(prhs[9]);
    cutoff = mxGetScalar(prhs[10]);
    seed = static_cast<unsigned int>(mxGetScalar(prhs[11]));

    openMat(matFile,para_file);

    para = matGetVariable(matFile,"para");
    assert(mxGetClassID(para) == mxSTRUCT_CLASS);
    getFieldVar(tau_ed, para, 0, "tau_e");
    getFieldVar(tau_er, para, 0, "tau_er");
    getFieldVar(tau_id, para, 0, "tau_i");
    getFieldVar(tau_ir, para, 0, "tau_ir");
    getFieldVar(vNa, para, 0, "vNa");
    getFieldVar(vK, para, 0, "vK");
    getFieldVar(vE, para, 0, "vE");
    getFieldVar(vI, para, 0, "vI");
    
    getIthElementOfFieldArray(gNa,para,ith,"gNa");
    getIthElementOfFieldArray(gK,para,ith,"gK");
    getIthElementOfFieldArray(gLeak,para,ith,"gLeak");
    getIthElementOfFieldArray(vLeak,para,ith,"vLeak");
    getIthElementOfFieldArray(vT,para,ith,"vT");
    getIthElementOfFieldArray(vRest,para,ith,"vRest");
    getIthElementOfFieldArray(DeltaT,para,ith,"DeltaT");
    getIthElementOfFieldArray(S,para,ith,"S");

    mxDestroyArray(para);
    closeMat(matFile,para_file);
    double pairs[3*2+2+3] = {gNa,vNa,gK,vK,gLeak,vLeak,vE,vI,vT,vRest,DeltaT};
                          //  0,  1, 2, 3,    4,    5, 6, 7, 8,   9,    10    
    for(i=0;i<neuroLib.nE;i++) {
        neuroLib.fE[i] = neuroLib.fE[i]/S;
        fStrength.push_back(neuroLib.fE[i]);
    }
    for(i=0;i<neuroLib.nI;i++) {
        neuroLib.fI[i] = neuroLib.fI[i]/S;
        fStrength.push_back(-neuroLib.fI[i]);
    }
    neuron.initialize(fStrength,neuroLib.nE,neuroLib.nI,0,seed,0,1);
    
    double tstep = neuroLib.tstep;
    nt = static_cast<unsigned int>(run_t/tstep)+1;
    simV.assign(nt,vinit);
    biV.assign(nt,vinit);
    liV.assign(nt,vinit);
    if (model==0) {
        m.assign(nt,m_inf(vinit,vT));
        n.assign(nt,n_inf(vinit,vT));
        h.assign(nt,h_inf(vinit,vT));
    }
    gE.assign(nt,0);
    gI.assign(nt,0);
    // get external inputs
    rE = rE/1000;
    rI = rI/1000; 
    while (neuron.status) neuron.getNextInput(rE,rE,rI,rI,run_t);

    // sim
    clock_gettime(clk_id,&tpS);
    neuron.vReset = vRest;
    unsigned int nc = 0;
    switch (model)
    {   
        case 0: //HH
            neuron.vThres = vRest + (vT -vRest)*1.5;
            cout << "HH " << endl;
            if (rk==2)
                nc = RK2_HH(simV,m,n,h,gE,gI,neuron,pairs,tau_er,tau_ed,tau_ir,tau_id,nt,tstep,cutoff,tsp_sim);
            if (rk==4)
                nc = RK4_HH(simV,m,n,h,gE,gI,neuron,pairs,tau_er,tau_ed,tau_ir,tau_id,nt,tstep,cutoff,tsp_sim);
            break;
        case 1: //EIF
            cout << "EIF " << endl;
            neuron.vThres = vNa;
            if (rk==2)
                nc = RK2_IF(simV,gE,gI,neuron,pairs,tau_er,tau_ed,tau_ir,tau_id,nt,tstep,cutoff,1,tsp_sim);
            if (rk==4)
                nc = RK4_IF(simV,gE,gI,neuron,pairs,tau_er,tau_ed,tau_ir,tau_id,nt,tstep,cutoff,1,tsp_sim);
            break;
        case 2: //IF
            cout << "IF " << endl;
            neuron.vThres = vT;
            if (rk==2)
                nc = RK2_IF(simV,gE,gI,neuron,pairs,tau_er,tau_ed,tau_ir,tau_id,nt,tstep,cutoff,0,tsp_sim);
            if (rk==4)
                nc = RK4_IF(simV,gE,gI,neuron,pairs,tau_er,tau_ed,tau_ir,tau_id,nt,tstep,cutoff,0,tsp_sim);
            break;
    }
    if (rk!=2 && rk!=4) cout << "rk2 or 4 only" << endl;
    else cout << "spikes: " << nc << endl;
    clock_gettime(clk_id,&tpE);
    cpu_t_sim = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
    cout << "sim ended, took " << cpu_t_sim << "s" << endl;

    // bilinear
    clock_gettime(clk_id,&tpS);
    switch (model) {   
        case 0: //HH
            neuron.vThres = vRest + (vT -vRest)*1.5;
            bilinear_HH(biV, neuroLib, neuron, run_t, ignore_t, tsp_bi);
            break;
        case 1: //EIF
            neuron.vThres = vNa;
            //bilinear_EIF();
            break;
        case 2: //IF
            neuron.vThres = vT;
            //bilinear_IF();
            break;
    }
    clock_gettime(clk_id,&tpE);
    cpu_t_bilinear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
    cout << "bilinear est. ended, took " << cpu_t_bilinear << "s" << endl;

    // linear 
    clock_gettime(clk_id,&tpS);
    switch (model) {   
        case 0: //HH
            neuron.vThres = vRest + (vT -vRest)*1.5;
            linear_HH(liV, neuroLib, neuron, run_t, ignore_t, tsp_li);
            break;
        case 1: //EIF
            neuron.vThres = vNa;
            //bilinear_EIF();
            break;
        case 2: //IF
            neuron.vThres = vT;
            //bilinear_IF();
            break;
    } 
    clock_gettime(clk_id,&tpE);
    cpu_t_linear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
    cout << "linear est. ended, took " << cpu_t_linear << "s" << endl;
    
    plhs[0] = mxCreateDoubleScalar(static_cast<double>(cpu_t_sim));
    plhs[1] = mxCreateDoubleScalar(static_cast<double>(cpu_t_bilinear));

    cout << " copy traces to output matrix " << endl;
    double *outTin, *outInID, *output;
    vector<double>* outputArray[8] = {&simV,&biV,&liV,&gE,&gI,&m,&n,&h};
    switch (model) {
        case 0:
            plhs[2] = mxCreateDoubleMatrix(nt, 8, mxREAL);
            output = mxGetPr(plhs[2]);
            for (j=0;j<8;j++) {
                for (i=0;i<nt;i++) {
                    *(output+j*nt+i) = outputArray[j]->at(i); 
                }   
            }
            break;
        default:
            plhs[2] = mxCreateDoubleMatrix(nt, 5, mxREAL);
            output = mxGetPr(plhs[2]);
            for (j=0;j<4;j++) {
                for (i=0;i<nt;i++) {
                    *(output+j*nt+i) = outputArray[j]->at(i); 
                }
            }
    }
    plhs[3] = mxCreateDoubleMatrix(neuron.tin.size(), 1, mxREAL);

    cout << " copy inputs to output matrix " << endl;
    outTin = mxGetPr(plhs[3]);
    for (i=0;i<neuron.tin.size();i++)
        outTin[i] = neuron.tin[i];
    //plhs[4] = mxCreateNumericMatrix(neuron.tin.size(), 1, mxINT64_CLASS, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(neuron.tin.size(), 1, mxREAL);
    outInID = mxGetPr(plhs[4]);
    for (i=0;i<neuron.tin.size();i++) {
        outInID[i] = static_cast<double>(neuron.inID[i]);
 //       cout << outInID[i] << endl;
    }
//    mexCallMATLAB(0,NULL,1,&(plhs[4]),"disp");
    neuroLib.clearLib();
    cout << "size: " <<  sizeof(size) << " bytes" << endl;
    cout << "size_b: "<< sizeof(size_b) << " bytes" << endl;
}
