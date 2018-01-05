#ifndef SINGLE_RK2_H
#define SINGLE_RK2_H
#include "NeuronStruct.h"
#include "math.h"
#include <vector>
#include <iostream>
#include "channels.h"
#include "typedefs.h"

void getCond2HH(std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &tin, std::vector<size_b> &inID,double tau_er, double tau_ed, double tau_ir, double tau_id, std::vector<double> &fStrength, double tstep) {
    const double conE = 1./(tau_er-tau_ed);
    const double conI = 1./(tau_ir-tau_id);
    double f,t;
    unsigned int i,j, nt = gE.size();

    for (j=0; j<tin.size(); j++) {
        f = fStrength[-inID[j]];
        if (f > 0)
            for (i=ceil(tin[j]/tstep); i<nt; i++) {
                t = i*tstep;
                //etr = exp(-(t-tin[j])/tau_er);
                gE[i] = gE[i] + f*conE*(exp(-(t-tin[j])/tau_er)-exp(-(t-tin[j])/tau_ed));
            }
        else
            for (i=ceil(tin[j]/tstep); i<nt; i++) {
                t = i*tstep;
                //etr = exp(-(t-tin[j])/tau_ir);
                gI[i] = gI[i] - f*conI*(exp(-(t-tin[j])/tau_ir)-exp(-(t-tin[j])/tau_id));
            }
    }
}

void getCond2IF(std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &hE, std::vector<double> &hI, std::vector<std::vector<double>> &tinBin, std::vector<std::vector<double>> &fBin, std::vector<double> &tin, std::vector<size_b> &inID, double tau_er, double tau_ed, double tau_ir, double tau_id, std::vector<double> &fStrength, double tstep) {
    const double conE = 1./(tau_er-tau_ed);
    const double conI = 1./(tau_ir-tau_id);
    double f,t,etr;
    unsigned int i,j, nt = gE.size();

    for (j=0; j<tin.size(); j++) {
        f = fStrength[-inID[j]];
        if (f > 0)
            for (i=ceil(tin[j]/tstep); i<nt; i++) {
                t = i*tstep;
                etr = exp(-(t-tin[j])/tau_er);
                gE[i] = gE[i] + f*conE*(etr-exp(-(t-tin[j])/tau_ed));
                hE[i] = hE[i] + f/tau_er*etr;
                tinBin[i].push_back(tin[j]);
                fBin[i].push_back(f);
            }
        else
            for (i=ceil(tin[j]/tstep); i<nt; i++) {
                t = i*tstep;
                etr = exp(-(t-tin[j])/tau_ir);
                gI[i] = gI[i] - f*conI*(etr-exp(-(t-tin[j])/tau_id));
                hI[i] = hI[i] - f/tau_ir*etr;
                tinBin[i].push_back(tin[j]);
                fBin[i].push_back(f);
            }
    }
}

void getCond2t(double gE, double gI, double hE, double hI, double &gEr, double &gIr, std::vector<double> &tinBin, std::vector<double> &fBin, double tsp, double t, double tau_er, double  tau_ed, double  tau_ir, double tau_id) {
    const double conE = 1./(tau_er-tau_ed);
    const double conI = 1./(tau_ir-tau_id);
    double etr,etd,c,dt;
    int i;
    etd = exp(-tsp/tau_ed);
    etr = exp(-tsp/tau_er);
    c = tau_er*conE*(etr-etd);
    gEr = gE * etd + c * hE;

    etd = exp(-tsp/tau_id);
    etr = exp(-tsp/tau_ir);
    c = tau_ir*conI*(etr-etd);
    gIr = gI * etd + c * hI;

    for (i=0; i<tinBin.size(); i++) {
        dt = tsp - (tinBin[i] - t);
        if (dt <= 0)
            break;
        if (fBin[i]>0) {
            etd = exp(-dt/tau_ed); 
            etr = exp(-dt/tau_er);
            gEr = gEr + fBin[i]*conE*(etr-etd);
        }
        else {
            etd = exp(-dt/tau_id); 
            etr = exp(-dt/tau_ir);
            gIr = gIr - fBin[i]*conI*(etr-etd);
        }
    }
}

unsigned int RK2_IF(std::vector<double> &v, std::vector<double> &gE, std::vector<double> &gI, Neuron neuron, double *pairs, double tau_er, double tau_ed, double tau_ir, double tau_id, size nt, double tstep, bool cutoff, bool eif, std::vector<double> &tsp) {
    unsigned int i, spikeCount = 0;
    double fk[2], vThres, tstep2 = tstep/2;
    double gEr, gIr;
    std::vector<double> hE(nt,0), hI(nt,0);
    std::vector<std::vector<double>> tinBin(nt,std::vector<double>());
    std::vector<std::vector<double>> fBin(nt,std::vector<double>());
    if (eif) vThres = pairs[1];
    else vThres = pairs[8];
    getCond2IF(gE,gI,hE,hI,tinBin,fBin,neuron.tin,neuron.inID,tau_er,tau_ed,tau_ir,tau_id,neuron.fStrength,tstep);
    for(i=1; i<nt; i++) {
        fk[0] = fk_IF(v[i-1],gE[i-1],gI[i-1],pairs,eif);
        fk[1] = fk_IF(v[i-1] + fk[0]*tstep,gE[i],gI[i],pairs,eif);
        v[i] = v[i-1] + tstep2*(fk[0] + fk[1]);

        if (v[i]> vThres) {
            spikeCount = spikeCount + 1;
            tsp.push_back((vThres-v[i-1])/(v[i]-v[i-1])*tstep);
            getCond2t(gE[i-1],gI[i-1],hE[i-1],hI[i-1],gEr,gIr,tinBin[i],fBin[i],tsp.back(),(i-1)*tstep,tau_er,tau_ed,tau_ir,tau_id);
            fk[0] = fk_IF(pairs[9],gEr,gIr,pairs,eif);
            fk[1] = fk_IF(pairs[9]+fk[0]*(tstep-tsp.back()),gE[i],gI[i],pairs,eif);
            v[i] = pairs[9] + (tstep-tsp.back())/2*(fk[0]+fk[1]);
            if (v[i] > vThres) std::cout << "use smaller tstep" << std::endl;
            if (cutoff) {
                //tsp = (i-1) * tstep + tsp; 
                break;
            }
        }
    }
    return spikeCount;
}	

unsigned int RK2_HH(std::vector<double> &v, std::vector<double> &m, std::vector<double> &n,std::vector<double> &h, std::vector<double> &gE,  std::vector<double> &gI, Neuron neuron, double *pairs, double tau_er, double tau_ed, double tau_ir, double tau_id, size nt, double tstep, std::vector<double> &tsp) {
    unsigned int i, spikeCount = 0;
    double mk[2],nk[2],hk[2],fk[2];
    double vT = pairs[8];
    double m_next, n_next, h_next, v_next, tstep2 = tstep/2.0;
    bool spiking = 0;
    getCond2HH(gE,gI,neuron.tin,neuron.inID,tau_er,tau_ed,tau_ir,tau_id,neuron.fStrength,tstep);
    for(i=1; i<nt; i++) {
        fk[0] = fk_HH(v[i-1],m[i-1],n[i-1],h[i-1],gE[i-1],gI[i-1],pairs); 
        v_next = v[i-1] + fk[0]*tstep;
        mk[0] = gatingm(v[i-1],m[i-1],vT);      m_next = m[i-1]+mk[0]*tstep;
        nk[0] = gatingn(v[i-1],n[i-1],vT);      n_next = n[i-1]+nk[0]*tstep; 
        hk[0] = gatingh(v[i-1],h[i-1],vT);      h_next = h[i-1]+hk[0]*tstep;
            
        fk[1] = fk_HH(v_next,m_next,n_next,h_next, gE[i],gI[i],pairs);
        mk[1] = gatingm(v_next,m_next,vT);   
        nk[1] = gatingn(v_next,n_next,vT);  
        hk[1] = gatingh(v_next,h_next,vT); 

        v[i] = v[i-1] + tstep2*(fk[0] + fk[1]);
        m[i] = m[i-1] + tstep2*(mk[0] + mk[1]);
        n[i] = n[i-1] + tstep2*(nk[0] + nk[1]);
        h[i] = h[i-1] + tstep2*(hk[0] + hk[1]);

        if (v[i]> neuron.vThres + (neuron.vThres-neuron.vReset)*2 && !spiking) {
            spikeCount = spikeCount + 1;
            if (v[i]<=v[i-1]) {
                tsp.push_back((i-1)*tstep);
            }
            spiking = 1;
        } else if (v[i] < neuron.vThres) spiking = 0;
    }
    return spikeCount;
}
#endif
