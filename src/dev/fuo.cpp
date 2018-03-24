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
#include "nNeuroSt.h"
#include "singleRK4.h"
//#include "singleRK2.h"
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
    ofstream tIncome_file, raster_file, data_file, cpu_file;
    ifstream cfg_file;
    string lib_file, para_file;
    mxArray *para;
    int ic;
    bool type[5];
    bool jumpy;
    bool shift;
    double cpu_t_sim, cpu_t_bilinear, cpu_t_linear, cpu_t_jbilinear;
    double rdHH,rgHH;
    clockid_t clk_id = CLOCK_PROCESS_CPUTIME_ID;
    struct timespec tpS, tpE;
    unsigned int ith, nt;
    vector<double> rE, rI;
    double run_t, ignore_t, tau_ed, tau_er, tau_id, tau_ir;
    double gNa, vNa, gK, vK, gLeak, vLeak, vT, vE, vI, vRest, DeltaT, S;
    double gM, gL, gT, vCA, vX, tau_max;
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
    vector<double> p, q, r, s, u;
    vector<double> pl, ql, rl, sl, ul;
    vector<double> pb, qb, rb, sb, ub;
    vector<size> crossb, crossl;
    unsigned int seed;
    bool spikeShape;
    bool test;
    vector<size_b>testID;
    vector<double>testTime;
    bool kVStyle, linear0;
    int afterCrossBehavior;
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
        ("seed,s",po::value<unsigned int>(&seed)->default_value(static_cast<unsigned int>(std::time(NULL))),"seeding")
		("theme",po::value<string>(&theme),"parameter file")
		("ith,i", po::value<unsigned int>(&ith)->default_value(1), "i(>0) th neuron")
		("rE", po::value<vector<double>>()->multitoken()->composing(), "Exc poisson rate")
		("rI", po::value<vector<double>>()->multitoken()->composing(), "Inh poisson rate")
		("testTime", po::value<vector<double>>()->multitoken()->composing(), "test input time array")
		("testID", po::value<vector<size_b>>()->multitoken()->composing(), "test input Strength ID array")
		("run_t,t", po::value<double>(&run_t), "sim time")
        ("vinit,v", po::value<unsigned int>(&vinit), "vinit")
        ("tref",po::value<double>(&tref),"refractory period")
        ("rdHH",po::value<double>(&rdHH),"drop HH threshold")
        ("rgHH",po::value<double>(&rgHH),"go HH threshold")
        ("afterCrossBehavior",po::value<int>(&afterCrossBehavior)->default_value(2)," 0:skip, 1:linear, 2:bilinear")
        ("spikeShape",po::value<bool>(&spikeShape)->default_value(true)," if false, crossing is spiking")
        ("kVStyle",po::value<bool>(&kVStyle)->default_value(true)," if false, crossing is spiking")
        ("linear0",po::value<bool>(&linear0)->default_value(true)," if false, linear est. is voltage dependent")
        ("test",po::value<bool>(&test)," if true, use test input")
        ("jumpy",po::value<bool>(&jumpy)->default_value(true)," if true, use jumpy bilinear")
        ("shift",po::value<bool>(&shift)," if true,round up input time to nearest time step")
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
		return 1;
	}
	cfg_file.close();
   
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
   
    ith = ith-1; 
    getIthElementOfFieldArray(gNa,para,ith,"gNa");
    getIthElementOfFieldArray(gK,para,ith,"gK");
    getIthElementOfFieldArray(gLeak,para,ith,"gLeak");
    getIthElementOfFieldArray(vLeak,para,ith,"vLeak");
    getIthElementOfFieldArray(vT,para,ith,"vT");
    getIthElementOfFieldArray(vRest,para,ith,"vRest");
    getIthElementOfFieldArray(DeltaT,para,ith,"DeltaT");
    getIthElementOfFieldArray(gM,para,ith,"gM");
    getIthElementOfFieldArray(gL,para,ith,"gL");
    getIthElementOfFieldArray(gT,para,ith,"gT");
    getIthElementOfFieldArray(vCA,para,ith,"vCA");
    getIthElementOfFieldArray(vX,para,ith,"vX");
    getIthElementOfFieldArray(tau_max,para,ith,"tau_max");
    getIthElementOfFieldArray(ic,para,ith,"type");
    ic = ic - 1;
    getIthElementOfFieldArray(S,para,ith,"S");
    mxDestroyArray(para);

    mxArray* tmp;
    tmp = matGetVariable(matFile,"bool");
    cout << "channel bool: ";
    for (int i=0; i<5; i++) {
        type[i] = static_cast<bool>(*(mxGetPr(tmp)+i*5+ic));
        if (i<4){
            cout << type[i] << ", ";
        } else {
            cout << type[i] << endl;
        }
    }
    mxDestroyArray(tmp);
    closeMat(matFile,para_file.c_str());

    neuroLib.readLib(lib_file.c_str());
    
    double pairs[17] = {gNa,vNa,gK,vK,gLeak,vLeak,vE,vI,vT,vRest,DeltaT,gM,gL,gT,vCA,vX,tau_max};
    // pairs: gNa,vNa,gK,vK,gLeak,vLeak,vE,vI,
    //         0   1   2  3  4    5      6  7 
    //        vT,vRest,DeltaT,
    //         8   9    10
    //        gM,gL,gT,vCA,Vx,tau_max
    //        11 12 13  14 15   16
    theme = theme + "-" + to_string(seed);
	raster_file.open("Raster-" + theme + ".bin", ios::out|ios::binary);
	tIncome_file.open("tIn-" + theme + ".bin", ios::out|ios::binary);
	data_file.open("Data-" + theme + ".bin", ios::out|ios::binary);
	cpu_file.open("cpuTime-" + theme + ".bin", ios::out|ios::binary);
    if (!raster_file) {
        cout << " failed to open file for writing" << endl;
        return 1;
    }
    if (!tIncome_file) {
        cout << " failed to open file for writing" << endl;
        return 1;
    }
    if (!data_file) {
        cout << " failed to open file for writing" << endl;
        return 1;
    }
    if (!cpu_file) {
        cout << " failed to open file for writing" << endl;
        return 1;
    }
    for(size i=0;i<neuroLib.nE;i++) {
        neuroLib.fE[i] = neuroLib.fE[i]/S;
        fStrength.push_back(neuroLib.fE[i]);
    }
    for(size i=0;i<neuroLib.nI;i++) {
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
    double vCrossl = vRest + (vT -vRest)*rgHH;
    double vBackl = vRest + (vT -vRest)*rdHH;
    double vCrossb = vRest + (vT -vRest)*rgHH;
    double vBackb = vRest + (vT -vRest)*rdHH;
    cout << " linear -> HH  " << vCrossl << endl;
    cout << " HH -> linear " << vBackl << endl;
    size plchldr_size0,plchldr_size1;
    double plchldr_double;
    nt = static_cast<unsigned int>(run_t/tstep)+1;
    cout << tstep << " ms -> " << nt << "steps " <<endl;
    cout << "initial voltage " << neuroLib.vRange[vinit] << endl;  
    if (test) {
        rEl = 1;
        if (vm.count("testTime") && vm.count("testID")) {
            testTime = vm["testTime"].as<vector<double>>();
            testID = vm["testID"].as<vector<size_b>>();
        } else {
            cout << " no test input provided" << endl;
            return 0;
        }
    }
    for (int k=0; k<rEl; k++) {
        cout << "============== " << k+1 << " ==============" << endl;
        simV.assign(nt,neuroLib.vRange[vinit]);
        biV.assign(nt,neuroLib.vRange[vinit]);
        liV.assign(nt,neuroLib.vRange[vinit]);
        m.reserve(nt);
        n.reserve(nt);
        h.reserve(nt);
        p.reserve(nt);
        q.reserve(nt);
        r.reserve(nt);
        s.reserve(nt);
        u.reserve(nt);
        ml.reserve(nt);
        nl.reserve(nt);
        hl.reserve(nt);
        pl.reserve(nt);
        ql.reserve(nt);
        rl.reserve(nt);
        sl.reserve(nt);
        ul.reserve(nt);
        mb.reserve(nt);
        nb.reserve(nt);
        hb.reserve(nt);
        pb.reserve(nt);
        qb.reserve(nt);
        rb.reserve(nt);
        sb.reserve(nt);
        ub.reserve(nt);
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
        vector<vector<double>> tPoi(neuroLib.nE+neuroLib.nI,vector<double>());
        if (test) {
            for (int it=0; it<testTime.size(); it++) {
                neuron.tin.push_back(testTime[it]);
                neuron.inID.push_back(testID[it]);
                tPoi[neuron.inID.back()].push_back(neuron.tin.back());
            }
        } else {
            double rEt= rE[k]/1000;
            double rIt= rI[k]/1000; 
            while (neuron.status) {
                if (neuron.getNextInput(rEt,rEt,rIt,rIt,run_t,shift,tstep)>=0) {
                    tPoi[neuron.inID.back()].push_back(neuron.tin.back());
                }
            }
        }
        for (size i=0;i<neuroLib.nE+neuroLib.nI; i++) {
            cout << i << ": {";
            for (size j=0;j<tPoi[i].size();j++) {
                cout << tPoi[i][j];
                if (j < tPoi[i].size()-1) {
                    cout << ", ";
                }
            }
            cout << "}" << endl;
        }
        cout << endl;
        if (jumpy) {
            nNS nn(
        }
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
        p.push_back(p_inf(neuroLib.vRange[vinit]));
        q.push_back(q_inf(neuroLib.vRange[vinit]));
        r.push_back(r_inf(neuroLib.vRange[vinit]));
        s.push_back(s_inf(neuroLib.vRange[vinit],vX));
        u.push_back(u_inf(neuroLib.vRange[vinit],vX));
        unsigned int nc = 0;
        neuron.vThres = vRest + 2*(vT -vRest);
        cout << "HH start" << endl;
        plchldr_size1 = 0;
        nc = RK4_HH(simV,m,n,h,p,q,r,s,u,gE,gI,hE,hI,neuron,pairs,type,tau_er,tau_ed,tau_ir,tau_id,nt,tstep,tsp_sim,false,0,plchldr_size0,plchldr_size1,plchldr_double);
        clock_gettime(clk_id,&tpE);
        cpu_t_sim = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "sim ended, took " << cpu_t_sim << "s" << endl;
        cout << "spikes: " << nc << endl;

        // bilinear
        nc = 0;
        clock_gettime(clk_id,&tpS);
        cout << " bilinear start " << endl;
        neuron.vThres = vRest + 2*(vT -vRest);
        nc = bilinear_HH(biV, gEb, gIb, hEb, hIb, mb, nb, hb, pb, qb, rb, sb, ub, crossb, neuroLib, neuron, run_t, ignore_t, tsp_bi, pairs, type, tau_er, tau_ed, tau_ir, tau_id, vCrossb, vBackb , neuron.tref, afterCrossBehavior, spikeShape, kVStyle);
        clock_gettime(clk_id,&tpE);
        cpu_t_bilinear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "bilinear est. ended, took " << cpu_t_bilinear << "s" << endl;
        cout << "spikes: " << nc << endl;
        
        // linear 
        nc = 0;
        cout << " linear start " << endl;
        clock_gettime(clk_id,&tpS);
        nc = linear_HH(liV, gEl, gIl, hEl, hIl, ml, nl, hl, pl, ql, rl, sl, ul, crossl, neuroLib, neuron, run_t, ignore_t, tsp_li, pairs, type, tau_er, tau_ed, tau_ir, tau_id, vCrossl, vBackl , neuron.tref, afterCrossBehavior, spikeShape,linear0);
        clock_gettime(clk_id,&tpE);
        cpu_t_linear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
        cout << "linear est. ended, took " << cpu_t_linear << "s" << endl;
        cout << "spikes: " << nc << endl;

        if (jumpy) {
            nc = 0;
            clock_gettime(clk_id,&tpS);
            cout << " bilinear start " << endl;
            nc = jbilinear_HH(
            clock_gettime(clk_id,&tpE);
            cpu_t_jbilinear = static_cast<double>(tpE.tv_sec-tpS.tv_sec) + static_cast<double>(tpE.tv_nsec - tpS.tv_nsec)/1e9;
            cout << "bilinear est. ended, took " << cpu_t_jbilinear << "s" << endl;
            cout << "spikes: " << nc << endl;
        }
       
        // data output
        writeTimeID(neuron.tin, neuron.inID, tIncome_file);
        
        vector<vector<double>> outTsp;
        outTsp.push_back(tsp_sim);
        outTsp.push_back(tsp_bi);
        outTsp.push_back(tsp_li);
        for (size i = 0; i< 3; i++) {
            writeTime(outTsp[i], raster_file);
        }
        outTsp.clear();

        cpu_file.write((char*)&(cpu_t_sim),sizeof(double));
        cpu_file.write((char*)&(cpu_t_bilinear), sizeof(double));
        cpu_file.write((char*)&(cpu_t_linear), sizeof(double));

        vector<double>* output[8] = {&simV,&biV,&liV,&gE,&gI,&m,&n,&h};
        for (size i=0;i<8;i++) {
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
        p.clear();
        q.clear();
        r.clear();
        s.clear();
        u.clear();
        ml.clear();
        nl.clear();
        hl.clear();
        pl.clear();
        ql.clear();
        rl.clear();
        sl.clear();
        ul.clear();
        mb.clear();
        nb.clear();
        hb.clear();
        pb.clear();
        qb.clear();
        rb.clear();
        sb.clear();
        ub.clear();
        simV.clear();
        biV.clear();
        liV.clear();
        crossl.clear();
        crossb.clear();
    }
    if (data_file.is_open())        data_file.close();
    if (raster_file.is_open())      raster_file.close();
    if (tIncome_file.is_open())     tIncome_file.close();
    if (cpu_file.is_open())         cpu_file.close();
    cout << "size: " <<  sizeof(size) << " bytes" << endl;
    cout << "size_b: "<< sizeof(size_b) << " bytes" << endl;
    neuroLib.clearLib();
	return 0;
}
