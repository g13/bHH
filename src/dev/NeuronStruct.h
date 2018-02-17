#ifndef NEURON_H
#define NEURON_H
#include "NeuronLibrary.h"
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <fstream>
#include <array>

//using std::vector

typedef struct NeuronStruct{
    size postN;
    size preN;
    double vThres, vReset, tref;
    unsigned int libID, nfE, nfI;
    // negative for inh
    std::vector<double> preStrength;
    std::vector<double> fStrength;
    //
    std::vector<size> postID;
    std::vector<size> preID;
    // 
    std::vector<size> IDatPost;
    //
    std::vector<double> tin;
    std::vector<size_b> inID; 
    std::vector<double> tinC;
    std::vector<size_b> CinID; 
    std::vector<double> tsp;
    size nOut;
    //
	std::mt19937_64 poiGenE, poiGenI;
	std::minstd_rand ranGenE, ranGenI, randE, randI;
    double tPoiI, tPoiE;
    bool Iready, Eready, ei, status, varInput, extE, extI;
    std::uniform_real_distribution<double> uniform0_1 = std::uniform_real_distribution<double>(0.0,1.0);;
    std::uniform_int_distribution<size_b> uniformE, uniformI;

    void initialize(std::vector<double> fStrength0, unsigned int nfE0, unsigned int nfI0, unsigned int libID0, unsigned int seed, size preN0, bool var) {
        libID = libID0;
        tPoiE = 0;
        tPoiI = 0;
        Iready = false;
        Eready = false;
        status = true;
        nOut = 0;
        preN = preN0;
        nfE = nfE0;
        nfI = nfI0;
        uniformE = std::uniform_int_distribution<size_b>(0,nfE-1);
        uniformI = std::uniform_int_distribution<size_b>(nfE,nfI+nfE-1);;
        if (nfE > 0) extE = 1;
        else extE = 0;
        if (nfI > 0) extI = 1;
        else extI = 0;
            
        varInput = var;
        //preID.assign(preID0.begin(),preID0.begin()+postN);
        fStrength.assign(fStrength0.begin(),fStrength0.end());
        //preStrength.assign(preStrength0.begin(),fStrength0.begin()+preN);
		std::array<unsigned int, std::mt19937_64::state_size> seqContainer;

		ranGenE.seed(seed);
		std::generate_n(seqContainer.data(), seqContainer.size(), std::ref(ranGenE));
		std::seed_seq seqE(seqContainer.begin(), seqContainer.end());
		poiGenE.seed(seqE);

		ranGenI.seed(seed+1);
		std::generate_n(seqContainer.data(), seqContainer.size(), std::ref(ranGenI));
		std::seed_seq seqI(seqContainer.begin(), seqContainer.end());
		poiGenI.seed(seqI);

        randE.seed(seed+2);
        randI.seed(seed+3);
    }

    double getNextInput(double rateE, double max_rateE, double rateI, double max_rateI, double t, bool shift, double tstep) {
        if (!Eready && extE) {
            if (rateE != max_rateE) {
                do tPoiE = tPoiE - log(uniform0_1(poiGenE))/max_rateE;
                while (uniform0_1(ranGenE) > rateE/max_rateE);
            } else tPoiE = tPoiE - log(uniform0_1(poiGenE))/max_rateE;
            if (shift){
                tPoiE = round(tPoiE/tstep)*tstep;
            }
        }
        if (!Iready && extI) {
            if (rateI != max_rateI) {
                do tPoiI = tPoiI - log(uniform0_1(poiGenI))/max_rateI;
                while (uniform0_1(ranGenI) > rateI/max_rateI);
            } else tPoiI = tPoiI - log(uniform0_1(poiGenI))/max_rateI;
            if (shift){
                tPoiI = round(tPoiE/tstep)*tstep;
            }
        }
        if (t < tPoiE && t < tPoiI) {
            status = 0;
            return -1;
        }
        if (tPoiI < tPoiE && extI) {
            Iready = 0;
            Eready = 1;
            tin.push_back(tPoiI);
            if (varInput) inID.push_back(uniformI(randI));
            else inID.push_back(2);
            return tPoiI;
        } else {
            if (extE) {
                Eready = 0;
                Iready = 1;
                tin.push_back(tPoiE);
                if (varInput) inID.push_back(uniformE(randE));
                else inID.push_back(1);
                return tPoiE;
            } else {
                status = 0;
                return -1;
            }
        }
    }
    void writeAndUpdateIn(size i, std::ofstream &tIncome_file) {
        
        // write 
        tIncome_file.write((char*)&i, sizeof(i));
        if (i!=0) {
            tIncome_file.write((char*)&(tin[0]), i * sizeof(tin[0]));
            tIncome_file.write((char*)&(inID[0]), i * sizeof(inID[0]));
        }
        size keeping = tin.size()-i;
        if (keeping) {
            tin.assign(tin.end()-keeping, tin.end());
            inID.assign(inID.end()-keeping, inID.end());
        } else {
            tin.clear();
            inID.clear();
        }
    }
    void writeAndUpdateOut(size i, std::ofstream &raster_file) {
        raster_file.write((char*)&i, sizeof(i));
        if (i!=0) 
            raster_file.write((char*)&(tsp[0]), i*sizeof(tsp[0]));
        nOut = nOut + i; 
        size keeping = tsp.size()-i;
        if (keeping) 
            tsp.assign(tsp.end()-keeping, tsp.end());
        else tsp.clear(); 
    }

    void setNextInput(double t, size ID){
        tinC.push_back(t);
        CinID.push_back(ID);
    }
    void clear() {
        tin.clear();
        inID.clear();
        tinC.clear();
        CinID.clear();
        tsp.clear();
        Iready = false;
        Eready = false;
        status = true;
        tPoiE = 0;
        tPoiI = 0;
    }
} Neuron;

void writeTime(std::vector<double> &t, std::ofstream &raster_file) {
    int i = t.size(); 
    raster_file.write((char*)&i, sizeof(i));
    if (i!=0) {
        raster_file.write((char*)&(t[0]), i*sizeof(t[0]));
    }
}

int writeTimeID(std::vector<double> &t, std::vector<size_b> &ID, std::ofstream &raster_file) {
    if (t.size() != ID.size()) {
        std::cout << " ID and t need to have the same size" << std::endl;
        return 0;
    }
    writeTime(t, raster_file);
    if (ID.size()>0) {
        raster_file.write((char*)&(ID[0]), ID.size()*sizeof(ID[0]));
    }
}
#endif
