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
    double *rE, *rI, run_t, ignore_t, tau_ed, tau_er, tau_id, tau_ir, rEt, rIt;
    double gNa, vNa, gK, vK, gLeak, vLeak, vT, vinit, vE, vI, vRest, DeltaT, S;
    double tsp_sim, tsp_bi, tsp_li, tref;
    MATFile *matFile;
    NeuroLib neuroLib;
    Neuron neuron;
    vector<double> fE,fI,fStrength;
    vector<double> simV, gE, gI, m, n, h, biV, liV;
    vector<double> gEb, gEl, gIb, gIl, ml, nl, hl, mb, nb, hb;
    vector<double> hEb, hEl, hIb, hIl, hE, hI;
    vector<size> crossb, crossl;
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

    mwSize rEl = mxGetN(prhs[3]);
    mwSize rIl = mxGetN(prhs[4]);
    assert(rIl==rEl);
    rE = new double[rEl];
    memcpy((void *) rE, (void *)(mxGetPr(prhs[3])),rEl*sizeof(double));
    cout << "rE: ";
    for (i=0;i<rEl-1;i++) cout << rE[i] << ", ";
    cout << rE[rEl-1] << endl;
    rI = new double[rIl];
    memcpy((void *) rI, (void *)(mxGetPr(prhs[4])),rIl*sizeof(double));
    cout << "rI: ";
    for (i=0;i<rIl-1;i++) cout << rI[i] << ", ";
    cout << rI[rIl-1] << endl;

    model = static_cast<unsigned int>(mxGetScalar(prhs[5]));
    run_t = mxGetScalar(prhs[6]);
    ignore_t = mxGetScalar(prhs[7]);
    vinit = mxGetScalar(prhs[8]);
    rk = mxGetScalar(prhs[9]);
    cutoff = mxGetScalar(prhs[10]);
    seed = static_cast<unsigned int>(mxGetScalar(prhs[11]));
    tref = mxGetScalar(prhs[12]);
    

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
    double tstep = neuroLib.tstep;
    neuron.initialize(fStrength,neuroLib.nE,neuroLib.nI,0,seed,0,1);
    neuron.tref = tref/tstep;
    
    nt = static_cast<int>(run_t/tstep)+1;
    vector<double>* outputArray[8];
    double *outTin, *outInID, *output;
    
    plhs[0] = mxCreateDoubleMatrix(6, rEl, mxREAL);
    int fiveOrEight;
    if (model == 0)
        fiveOrEight = 8;
    else 
        fiveOrEight = 5;

    mwSize dims[3];
    dims[0] = nt;
    dims[1] = fiveOrEight;
    dims[2] = rEl;
    plhs[1] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    //tmp = new mxArray*[rEl];
    plhs[2] = mxCreateCellMatrix(rEl,1);
    plhs[3] = mxCreateCellMatrix(rEl,1);
    unsigned int nc;
    double rHH = 1.3;
    double rLinear = 0.7;
    double vCrossl = vRest + (vT -vRest)*rLinear;
    double vBackl = vRest + (vT -vRest)*rLinear;
    cout << " linear -> HH  " << vCrossl << endl;
    cout << " HH -> linear " << vBackl << endl;
    size plchldr_size0,plchldr_size1;
    double plchldr_double;
    cout << tstep << " ms -> " << nt << "steps " <<endl;
    cout << "initial voltage " << vinit << endl;  
    for (k=0;k<rEl;k++) {
        cout << "============== " << k+1 << " ==============" << endl;
        simV.assign(nt,vinit);
        biV.assign(nt,vinit);
        liV.assign(nt,vinit);
        if (model==0) {
            m.reserve(nt);
            n.reserve(nt);
            h.reserve(nt);
            ml.reserve(nt);
            nl.reserve(nt);
            hl.reserve(nt);
            mb.reserve(nt);
            nb.reserve(nt);
            hb.reserve(nt);
        }
        gE.reserve(nt);
        gI.reserve(nt);
        gEb.reserve(nt);
        gIb.reserve(nt);
        gEl.reserve(nt);
        gIl.reserve(nt);
        hE.reserve(nt);
        hI.reserve(nt);
        hEb.reserve(nt);
        hIb.reserve(nt);
        hEl.reserve(nt);
        hIl.reserve(nt);
        // get external inputs
        rEt = rE[k]/1000;
        rIt = rI[k]/1000; 
        while (neuron.status) neuron.getNextInput(rEt,rEt,rIt,rIt,run_t);

        // sim
        clock_gettime(clk_id,&tpS);
        neuron.vReset = vRest;
        gE.push_back(0);
        gI.push_back(0);
        hE.push_back(0);
        hI.push_back(0);
        m.push_back(m_inf(vinit,vT));
        n.push_back(n_inf(vinit,vT));
        h.push_back(h_inf(vinit,vT));
        nc = 0;
        switch (model)
        {   
            case 0: //HH
                neuron.vThres = vRest + (vT -vRest)*rHH;
                cout << "HH " << endl;
                if (rk==2)
                    nc = RK2_HH(simV,m,n,h,gE,gI,neuron,pairs,tau_er,tau_ed,tau_ir,tau_id,nt,tstep,tsp_sim);
                if (rk==4) {
                    nc = RK4_HH(simV,m,n,h,gE,gI,hE,hI,neuron,pairs,tau_er,tau_ed,tau_ir,tau_id,nt,tstep,tsp_sim,false,0,plchldr_size0,plchldr_size1,plchldr_double);
                    assert(plchldr_size1 == neuron.tin.size()-1);
                }
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
        clock_gettime(clk_id,&tpE);
        cpu_t_sim = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "sim ended, took " << cpu_t_sim << "s" << endl;
        if (rk!=2 && rk!=4) cout << "rk2 or 4 only" << endl;
        else cout << "spikes: " << nc << endl;

        // bilinear
        clock_gettime(clk_id,&tpS);
        nc = 0;
        /*switch (model) {   
            case 0: //HH
                neuron.vThres = vRest + (vT -vRest)*rHH;
                nc = bilinear_HH(biV, neuroLib, neuron, run_t, ignore_t, tsp_bi);
                break;
            case 1: //EIF
                neuron.vThres = vNa;
                //bilinear_EIF();
                break;
            case 2: //IF
                neuron.vThres = vT;
                //bilinear_IF();
                break;
        }*/
        clock_gettime(clk_id,&tpE);
        cpu_t_bilinear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "bilinear est. ended, took " << cpu_t_bilinear << "s" << endl;
        cout << "spikes: " << nc << endl;

        // linear 
        clock_gettime(clk_id,&tpS);
        nc = 0;
        switch (model) {   
            case 0: //HH
                linear_HH(liV, gEl, gIl, hEl, hIl, ml, nl, hl, crossl, neuroLib, neuron, run_t, ignore_t, tsp_li, pairs, tau_er, tau_ed, tau_ir, tau_id, vCrossl, vBackl , neuron.tref, rk);
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
        cout << "spikes: " << nc << endl;

        cout << " copy traces to output matrix " << endl;
        output = mxGetPr(plhs[0]);
        *(output + k*6 + 0) = cpu_t_sim;
        *(output + k*6 + 1) = cpu_t_bilinear;
        *(output + k*6 + 2) = cpu_t_linear;
        *(output + k*6 + 3) = tsp_sim;
        *(output + k*6 + 4) = tsp_bi;
        *(output + k*6 + 5) = tsp_li;

        outputArray[0] = &(simV);
        outputArray[1] = &(biV);
        outputArray[2] = &(liV);
        outputArray[3] = &(gE);
        outputArray[4] = &(gI);
        outputArray[5] = &(m);
        outputArray[6] = &(n);
        outputArray[7] = &(h);

        output = mxGetPr(plhs[1]);
        for (j=0;j<fiveOrEight;j++) {
            for (i=0;i<nt;i++) {
                *(output+k*fiveOrEight*nt+j*nt+i) = outputArray[j]->at(i); 
            }   
        }
        cout << " copy inputs to output matrix " << endl;
        tmp = mxCreateDoubleMatrix(neuron.tin.size(),1,mxREAL);
        outTin = mxGetPr(tmp);
        for (i=0;i<neuron.tin.size();i++) {
            *(outTin+i) = neuron.tin[i];
        }
        mxSetCell(plhs[2], (mwIndex) k, mxDuplicateArray(tmp));
        mxDestroyArray(tmp);

        tmp = mxCreateDoubleMatrix(neuron.tin.size(),1,mxREAL);
        outInID = mxGetPr(tmp);
        for (i=0;i<neuron.tin.size();i++) {
            *(outInID+i) = static_cast<double>(neuron.inID[i]);
        }
        mxSetCell(plhs[3], (mwIndex) k, mxDuplicateArray(tmp));
        mxDestroyArray(tmp);
        neuron.clear();
    
        gE.clear();
        gI.clear();
        gEb.clear();
        gIb.clear();
        gEl.clear();
        gIl.clear();
        hE.clear();
        hI.clear();
        hEb.clear();
        hIb.clear();
        hEl.clear();
        hIl.clear();
        m.clear();
        n.clear();
        h.clear();
        ml.clear();
        nl.clear();
        hl.clear();
        mb.clear();
        nb.clear();
    }
    neuroLib.clearLib();
    cout << "size: " <<  sizeof(size) << " bytes" << endl;
    cout << "size_b: "<< sizeof(size_b) << " bytes" << endl;
    delete []rE;
    delete []rI;
    mxFree(lib_file);
    cout << "freed" << endl;
    mxFree(para_file);
    cout << "freed" << endl;
}
