//#include <fenv.h>
#include <boost/program_options.hpp>
#include "mex.h"
//#include "boost_program_options_overload.h"
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
    ofstream tIncome_file, raster_file, data_file;
    ifstream cfg_file;
    string lib_file, para_file;
    mxArray *para ;
    double cpu_t_sim, cpu_t_bilinear, cpu_t_linear;
    double rLinear;
    clockid_t clk_id = CLOCK_PROCESS_CPUTIME_ID;
    struct timespec tpS, tpE;
    unsigned int ith,i,j,nt;
    vector<double> rE, rI;
    double run_t, ignore_t, tau_ed, tau_er, tau_id, tau_ir;
    double gNa, vNa, gK, vK, gLeak, vLeak, vT, vE, vI, vRest, DeltaT,S;
    vector<double> tsp_sim, tsp_bi, tsp_li;
    double tref;
    unsigned int vinit;
    MATFile *matFile;
    NeuroLib neuroLib;
    Neuron neuron;
    po::variables_map vm;
    vector<double> fE,fI,fStrength;
    vector<double> simV, gE, gI, m, n, h, biV, liV;
    vector<double> gEb, gEl, gIb, gIl, ml, nl, hl, mb, nb, hb;
    vector<double> hEb, hEl, hIb, hIl, hE, hI;
    vector<size> crossb, crossl;
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
		("rE", po::value<vector<double>>()->multitoken()->composing(), "Exc poisson rate")
		("rI", po::value<vector<double>>()->multitoken()->composing(), "Inh poisson rate")
		("run_t,t", po::value<double>(&run_t), "sim time")
        ("vinit,v", po::value<unsigned int>(&vinit), "sim time")
        ("tref",po::value<double>(&tref),"refractory period")
        ("rLinear",po::value<double>(&rLinear),"linear->HH threshold")
		("ignore_t", po::value<double>(&ignore_t), "ingore time while applying bilinear rules");

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
    rE = vm["rE"].as<vector<double>>();
    rI = vm["rI"].as<vector<double>>();
    int rEl = rE.size();
    int rIl = rI.size();
    double rHH = 1.0;
    double rLinear = 0.8;
    double rBilinear = rLinear;
    double vCrossl = vRest + (vT -vRest);
    double vBackl = vRest + (vT -vRest)*rLinear;
    double vCrossb = vRest + (vT -vRest);
    double vBackb = vRest + (vT -vRest)*rBilinear;
    cout << " linear -> HH  " << vCrossl << endl;
    cout << " HH -> linear " << vBackl << endl;
    size plchldr_size0,plchldr_size1;
    double plchldr_double;
    nt = static_cast<unsigned int>(run_t/tstep)+1;
    cout << tstep << " ms -> " << nt << "steps " <<endl;
    cout << "initial voltage " << neuroLib.vRange[vinit] << endl;  
    for (int k=0; k<rEl; k++) {
        cout << "============== " << k+1 << " ==============" << endl;
        simV.assign(nt,neuroLib.vRange[vinit]);
        biV.assign(nt,neuroLib.vRange[vinit]);
        liV.assign(nt,neuroLib.vRange[vinit]);
        m.reserve(nt);
        n.reserve(nt);
        h.reserve(nt);
        ml.reserve(nt);
        nl.reserve(nt);
        hl.reserve(nt);
        mb.reserve(nt);
        nb.reserve(nt);
        hb.reserve(nt);
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
        double rEt= rE[k]/1000;
        double rIt= rI[k]/1000; 
        if (rE == 0.0) neuron.extE = false;
        if (rI == 0.0) neuron.extI = false;
        while (neuron.status) neuron.getNextInput(rEt,rEt,rIt,rIt,run_t);

        // HH sim 
        clock_gettime(clk_id,&tpS);
        neuron.vReset = vRest;
        gE.push_back(0);
        gI.push_back(0);
        hE.push_back(0);
        hI.push_back(0);
        m.push_back(m_inf(neuroLib.vRange[vinit],vT));
        n.push_back(n_inf(neuroLib.vRange[vinit],vT));
        h.push_back(h_inf(neuroLib.vRange[vinit],vT));
        unsigned int nc = 0;
        neuron.vThres = vRest + 2*(vT -vRest)*rHH;
        cout << "HH start" << endl;
        nc = RK4_HH(simV,m,n,h,gE,gI,hE,hI,neuron,pairs,tau_er,tau_ed,tau_ir,tau_id,nt,tstep,tsp_sim,false,0,plchldr_size0,plchldr_size1,plchldr_double);
        clock_gettime(clk_id,&tpE);
        cpu_t_sim = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "sim ended, took " << cpu_t_sim << "s" << endl;
        cout << "spikes: " << nc << endl;

        // bilinear
        nc = 0;
        clock_gettime(clk_id,&tpS);
        cout << " bilinear start " << endl;
        neuron.vThres = vRest + (vT -vRest)*rHH;
        nc = bilinear_HH(biV, gEb, gIb, hEb, hIb, mb, nb, hb, crossb, neuroLib, neuron, run_t, ignore_t, tsp_bi, pairs, tau_er, tau_ed, tau_ir, tau_id, vCrossb, vBackb , neuron.tref, afterCrossBehavior, spikeShape, kVStyle);
        clock_gettime(clk_id,&tpE);
        cpu_t_bilinear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "bilinear est. ended, took " << cpu_t_bilinear << "s" << endl;
        cout << "spikes: " << nc << endl;
        
        // linear 
        nc = 0;
        cout << " linear start " << endl;
        clock_gettime(clk_id,&tpS);
        nc = linear_HH(liV, gEl, gIl, hEl, hIl, ml, nl, hl, crossl, neuroLib, neuron, run_t, ignore_t, tsp_li, pairs, tau_er, tau_ed, tau_ir, tau_id, vCrossl, vBackl , neuron.tref, afterCrossBehavior, spikeShape);
        clock_gettime(clk_id,&tpE);
        cpu_t_linear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "linear est. ended, took " << cpu_t_linear << "s" << endl;
        cout << "spikes: " << nc << endl;
       
        // data output
        writeTimeID(neuron.tin, neuron.inID tIncome_file);
        
        vector<double> outTsp;
        outTsp.push_back(tsp_sim);
        outTsp.push_back(tsp_bi);
        outTsp.push_back(tsp_li);
        for (i = 0; i< 3; i++) {
            writeTime(outTsp[i], raster_file);
        }
        outTsp.clear();

        vector<double>* output[8] = {&simV,&biV,&liV,&gE,&gI,&m,&n,&h};
        for (i=0;i<8;i++) {
            data_file.write((char*)&(output[i]->at(0)), nt*sizeof(double));
        }
        neuron.clear();
        tsp_sim.clear();
        tsp_bi.clear();
        tsp_li.clear();
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
        simV.clear();
        biV.clear();
        liV.clear();
    }
    neuroLib.clearLib();
    cout << "size: " <<  sizeof(size) << " bytes" << endl;
    cout << "size_b: "<< sizeof(size_b) << " bytes" << endl;
    if (data_file.is_open())        data_file.close();
    if (raster_file.is_open())      raster_file.close();
    if (tIncome_file.is_open())     tIncome_file.close();
}
