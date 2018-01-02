//#include <fenv.h>
#include "mex.h"
#include <boost/program_options.hpp>
#include <cstring>
//#include "boost_program_options_overload.h"
#include "matFunc.h"
#include "time.h"
#include "channels.h"
#include "NeuronLibrary.h"
#include "NeuronStruct.h"
#include "singleRK4.h"
#include "singleRK2.h"
#include "singleBilinear.h"
#include "typedefs.h"
#include <fstream>
#include <iostream>

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ofstream;
using std::ifstream;
using std::to_string;
using std::ios;
namespace po = boost::program_options;

int main(int argc, char **argv)
{
    //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    // prhs[] ~ lib_file, para_file, ith, rE, rI, model, run_t, ignore_t, vinit
    ofstream tIncome_file, raster_file, data_file;
    ifstream cfg_file;
    string lib_file, para_file;
    mxArray *para ;
    double cpu_t_sim, cpu_t_bilinear, cpu_t_linear;
    double rLinear;
    clockid_t clk_id = CLOCK_PROCESS_CPUTIME_ID;
    struct timespec tpS, tpE;
    unsigned int ith,i,j,k,model,nt,rk;
    double rE, rI, run_t, ignore_t, tau_ed, tau_er, tau_id, tau_ir;
    double gNa, vNa, gK, vK, gLeak, vLeak, vT, vE, vI, vRest, DeltaT,S;
    double tsp_sim, tsp_simq, tsp_bi, tsp_li, tref;
    double tmpsum;
    unsigned int vinit;
    MATFile *matFile;
    NeuroLib neuroLib;
    Neuron neuron;
    po::variables_map vm;
    vector<double> fE,fI,fStrength;
    //vector<double> simVq, gEq, gIq, hEq, hIq, mq, nq, hq;
    vector<double> simV, gE, gI, hE, hI, m, n, h;
    vector<double> biV, gEb, gIb, hEb, hIb, mb, nb, hb;
    vector<double> liV;
    vector<double> gEl;
    vector<double> gIl;
    vector<double> hEl;
    vector<double> hIl;
    vector<double> ml;
    vector<double> nl;
    vector<double> hl;
    vector<size> crossb, crossl;
    bool cutoff;
	unsigned int seed0 = static_cast<unsigned int>(std::time(NULL));
    unsigned int seed;
    string config_file;
    string theme;

	po::options_description cmd_line_options("cmd line options"),
        config_file_options("config file"),
        Generic("options");
    cmd_line_options.add_options()
        ("config_file,c", po::value<string>(&config_file)->default_value("test.cfg"),"config filename");
	Generic.add_options()
		("lib_file,L",po::value<string>(&lib_file),"neuron library")
		("para_file,p",po::value<string>(&para_file),"parameter file")
        ("seed,s",po::value<unsigned int>(&seed),"seeding")
		("theme",po::value<string>(&theme),"parameter file")
		("ith,i", po::value<unsigned int>(&ith)->default_value(1), "i th neuron")
		("rE", po::value<double>(&rE), "Exc poisson rate")
		("rI", po::value<double>(&rI), "Exc poisson rate")
		("model,m", po::value<unsigned int>(&model), "model 0-HH 1-EIF 2-IF")
		("run_t,t", po::value<double>(&run_t), "sim time")
        ("vinit,v", po::value<unsigned int>(&vinit), "sim time")
        ("rk,r",po::value<unsigned int>(&rk), "order of rk method")
        ("cutoff",po::value<bool>(&cutoff), "terminate when spike")
        ("tref",po::value<double>(&tref),"refractory period")
        ("rLinear",po::value<double>(&rLinear),"linear->HH threshold")
		("ignore_t", po::value<double>(&ignore_t), "ingore time while applying bilinear rules");
        //("v",po::value<vector<double> >()->multitoken(),"vector test");

	cmd_line_options.add(Generic);
	config_file_options.add(Generic);
	store(po::parse_command_line(argc, argv, cmd_line_options), vm);
	notify(vm);
	cfg_file.open(config_file);
	if (cfg_file) {
		store(po::parse_config_file(cfg_file, config_file_options, true), vm);
		notify(vm);
	} else {
		cout << "cannot open the configuration file: " << config_file << endl;
		return 0;
	}
	cfg_file.close();
   
    neuroLib.readLib(lib_file.c_str());
    
    openMat(matFile,para_file.c_str());

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
    closeMat(matFile,para_file.c_str());
    double pairs[3*2+2+3] = {gNa,vNa,gK,vK,gLeak,vLeak,vE,vI,vT,vRest,DeltaT};
                          //  0,  1, 2, 3,    4,    5, 6, 7, 8,   9,    10    
    
    string prefix = to_string(static_cast<int>(rE))+"-"+to_string(static_cast<int>(rI))+"-";
    theme = "-" + theme;
	raster_file.open(prefix + "Raster" + theme + ".bin", ios::out|ios::binary);
	tIncome_file.open(prefix + "tIn" + theme + ".bin", ios::out|ios::binary);
	data_file.open(prefix + "Data" + theme + ".bin", ios::out|ios::binary);
    if (!seed) seed = seed0;
    for(i=0;i<neuroLib.nE;i++) {
        neuroLib.fE[i] = neuroLib.fE[i]/S;
        fStrength.push_back(neuroLib.fE[i]);
    }
    for(i=0;i<neuroLib.nI;i++) {
        neuroLib.fI[i] = neuroLib.fI[i]/S;
        fStrength.push_back(-neuroLib.fI[i]);
    }
    double tstep = neuroLib.tstep;
    cout << "using tstep: " << tstep << " ms" << endl;
    neuron.initialize(fStrength,neuroLib.nE,neuroLib.nI,0,seed,0,1);
    neuron.tref = tref;

    nt = static_cast<unsigned int>(run_t/tstep)+1;
    simV.assign(nt,neuroLib.vRange[vinit]);
    biV.assign(nt,neuroLib.vRange[vinit]);
    liV.assign(nt,neuroLib.vRange[vinit]);
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
    rE = rE/1000;
    rI = rI/1000; 
    if (rE == 0.0) neuron.extE = false;
    if (rI == 0.0) neuron.extI = false;
    while (neuron.status) neuron.getNextInput(rE,rE,rI,rI,run_t);

    // sim
    clock_gettime(clk_id,&tpS);
    neuron.vReset = vRest;
    gE.push_back(0);
    gI.push_back(0);
    hE.push_back(0);
    hI.push_back(0);
    m.push_back(m_inf(neuroLib.vRange[vinit],vT));
    n.push_back(n_inf(neuroLib.vRange[vinit],vT));
    h.push_back(h_inf(neuroLib.vRange[vinit],vT));
    double rHH = 1.3;
    double vCrossl = vRest + (vT -vRest)*rLinear;
    double vBackl = vRest + (vT -vRest)*rLinear - 1;
    cout << " linear -> HH  " << vCrossl << endl;
    cout << " HH -> linear " << vBackl << endl;
    size plchldr_size0,plchldr_size1;
    double plchldr_double;
    unsigned int nc;
    switch (model)
    {   
        case 0: //HH
            neuron.vThres = vRest + (vT -vRest)*rHH;
            cout << "HH " << endl;
            if (rk==2)
                nc = RK2_HH(simV,m,n,h,gE,gI,neuron,pairs,tau_er,tau_ed,tau_ir,tau_id,nt,tstep,tsp_sim);
            if (rk==4) {
                nc = RK4_HH(simV,m,n,h,gE,gI,hE,hI,neuron,pairs,tau_er,tau_ed,tau_ir,tau_id,nt,tstep,tsp_sim,false,0,plchldr_size0,plchldr_size1,plchldr_double);
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
    else {
        cout << "spikes: " << nc << endl;
        //cout << "spikes: " << nc1 << endl;
    }

    // bilinear
    clock_gettime(clk_id,&tpS);
    nc = 0;
    switch (model) {   
        case 0: //HH
            neuron.vThres = vRest + (vT -vRest)*1.5;
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
    }
    clock_gettime(clk_id,&tpE);
    cpu_t_bilinear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
    cout << "bilinear est. ended, took " << cpu_t_bilinear << "s" << endl;
    cout << "spikes: " << nc << endl;
    
    // linear 
    clock_gettime(clk_id,&tpS);
    nc = 0;
    switch (model) {   
        case 0: //HH
            nc = linear_HH(liV, gEl, gIl, hEl, hIl, ml, nl, hl, crossl, neuroLib, neuron, run_t, ignore_t, tsp_li, pairs, tau_er, tau_ed, tau_ir, tau_id, vCrossl, vBackl, neuron.tref, rk);
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
    

    vector<double>* output[8] = {&simV,&biV,&liV,&gE,&gI,&m,&n,&h};
    neuron.writeAndUpdateIn(neuron.tin.size(), tIncome_file);
    neuron.writeAndUpdateOut(neuron.tsp.size(), raster_file);
    switch (model) {
        case 0:
            for (i=0;i<8;i++)
                data_file.write((char*)&(output[i]->at(0)), nt*sizeof(double));
            break;
        default:
            for (i=0;i<5;i++)
                data_file.write((char*)&(output[i]->at(0)), nt*sizeof(double));
    }
    if (data_file.is_open())        data_file.close();
    if (raster_file.is_open())      raster_file.close();
    if (tIncome_file.is_open())     tIncome_file.close();
    neuroLib.clearLib();
    cout << "size: " <<  sizeof(size) << " bytes" << endl;
    cout << "size_b: "<< sizeof(size_b) << " bytes" << endl;
}
