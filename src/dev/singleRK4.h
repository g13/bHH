#ifndef SINGLE_RK4_H
#define SINGLE_RK4_H

#include "NeuronStruct.h"
#include "math.h"
#include <vector>
#include <iostream>
#include "channels.h"
#include "typedefs.h"

void getCond4t(double gE, double gI, double hE, double hI, double &gEr, double &gIr, double &gErmid, double &gIrmid, std::vector<double> &tinBin, std::vector<double> &fBin, double tsp, double t, double tau_er, double  tau_ed, double  tau_ir, double tau_id, double tstep) {
    const double conE = 1./(tau_er-tau_ed);
    const double conI = 1./(tau_ir-tau_id);
    double etr,etd,c,dt;
    int i;
    etd = exp(-tsp/tau_ed);
    etr = exp(-tsp/tau_er);
    c = tau_er*conE*(etr-etd);
    gEr = gE * etd + c * hE;
    hE = hE * etr;

    etd = exp(-tsp/tau_id);
    etr = exp(-tsp/tau_ir);
    c = tau_ir*conI*(etr-etd);
    gIr = gI * etd + c * hI;
    hI = hI * etr;

    for (i=0; i<tinBin.size(); i++) {
        dt = tsp - (tinBin[i] - t);
        if (dt <= 0)
            break;
        if (fBin[i]>0) {
            etd = exp(-dt/tau_ed); 
            etr = exp(-dt/tau_er);
            gEr = gEr + fBin[i]*conE*(etr-etd);
            hE = hE + fBin[i]/tau_er*etr;
        }
        else {
            etd = exp(-dt/tau_id); 
            etr = exp(-dt/tau_ir);
            gIr = gIr - fBin[i]*conI*(etr-etd);
            hI = hI - fBin[i]/tau_ir*etr;
        }
    }
    
    tsp = tsp + (tstep-tsp)/2;

    etd = exp(-tsp/tau_ed);
    etr = exp(-tsp/tau_er);
    c = tau_er*conE*(etr-etd);
    gErmid = gEr * etd + c * hE;

    etd = exp(-tsp/tau_id);
    etr = exp(-tsp/tau_ir);
    c = tau_ir*conI*(etr-etd);
    gIrmid = gIr * etd + c * hI;

    for (i=i; i<tinBin.size(); i++) {
        dt = tsp - (tinBin[i] - t);
        if (dt <= 0)
            break;
        if (fBin[i]>0) {
            etd = exp(-dt/tau_ed); 
            etr = exp(-dt/tau_er);
            gErmid = gErmid + fBin[i]*conE*(etr-etd);
        }
        else {
            etd = exp(-dt/tau_id); 
            etr = exp(-dt/tau_ir);
            gIrmid = gIrmid - fBin[i]*conI*(etr-etd);
        }
    }
}

void getCond4IF(std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &gEmid, std::vector<double> &gImid, std::vector<double> &hE, std::vector<double> &hI, std::vector<std::vector<double>> &tinBin, std::vector<std::vector<double>> &fBin, std::vector<double> &tin, std::vector<size_b> &inID, double tau_er, double tau_ed, double tau_ir, double tau_id, std::vector<double> &fStrength, double tstep, size nt) {
    const double conE = 1./(tau_er-tau_ed);
    const double conI = 1./(tau_ir-tau_id);
    double f,t,etr,tstep2 = tstep/2;
    unsigned int i,j;

    for (j=0; j<tin.size(); j++) {
        f = fStrength[inID[j]];
        if (f > 0)
            for (i=ceil(tin[j]/tstep); i<nt; i++) {
                t = i*tstep;
                etr = exp(-(t-tin[j])/tau_er);
                gE[i] = gE[i] + f*conE*(etr-exp(-(t-tin[j])/tau_ed));
                hE[i] = hE[i] + f/tau_er*etr;
                t = t + tstep2;
                gEmid[i] = gEmid[i] + f*conE*(exp(-(t-tin[j])/tau_er)-exp(-(t-tin[j])/tau_ed));
                tinBin[i].push_back(tin[j]);
                fBin[i].push_back(f);
            }
        else
            for (i=ceil(tin[j]/tstep); i<nt; i++) {
                t = i*tstep;
                etr = exp(-(t-tin[j])/tau_ir);
                gI[i] = gI[i] - f*conI*(etr-exp(-(t-tin[j])/tau_id));
                hI[i] = hI[i] - f/tau_ir*etr;
                t = t + tstep2;
                gImid[i] = gImid[i] - f*conI*(exp(-(t-tin[j])/tau_ir)-exp(-(t-tin[j])/tau_id));
                tinBin[i].push_back(tin[j]);
                fBin[i].push_back(f);
            }
    }
}

inline double alpha(double dt, double tau_d, double tau_r, double g, double &h) 
{ // note that h is a reference variable here, so the original variable will also change when the h in this function is modified 
	double  etd = exp(-dt/tau_d), etr = exp(-dt/tau_r);
	double c = tau_r/(tau_d-tau_r) * (etd-etr);
	g = g * etd + c * h;
	h = h * etr;
    return g;
}

void getCond4HH(std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &hE, std::vector<double> &hI, double &gEmid, double &gImid, std::vector<double> &tin, std::vector<size_b> &inID, size &is, size ns, double tau_er, double tau_ed, double tau_ir, double tau_id, std::vector<double> &fStrength, size i, double tstep, double t, bool extE, bool extI){
    double dt, dtE, dtI;
    double tstep2 = tstep/2;
    double f;
    bool noMid;
    unsigned int j;

    if (!extE) {
        gEmid = 0;
        gE.push_back(0);
        hE.push_back(0);
    } else {
        gE.push_back(gE[i-1]);
        hE.push_back(hE[i-1]);
        dtE = 0;
    }

    if (!extI) {
        gImid = 0;
        gI.push_back(0);
        hI.push_back(0);
    } else {
        gI.push_back(gI[i-1]);
        hI.push_back(hI[i-1]);
        dtI = 0;
    }
            
    if (!extE && !extI) return;
    noMid = true;
    if (is < ns) {
        while (tin[is] < t + tstep) {
            dt = tin[is] - t;
            if (dt > tstep2 && noMid) {
                if (extE) {
                    gEmid = alpha(tstep2-dtE,tau_ed,tau_er,gE[i],hE[i]);
                    gE[i] = gEmid;
                    dtE = tstep2;
                }
                if (extI) {
                    gImid = alpha(tstep2-dtI,tau_id,tau_ir,gI[i],hI[i]);
                    gI[i] = gImid;
                    dtI = tstep2;
                }
                noMid = false;
            }
             
            f = fStrength[inID[is]];
            if (f>0) {
                gE[i] = alpha(dt - dtE,tau_ed,tau_er,gE[i],hE[i]);
                dtE = dt;
                hE[i] = hE[i] + f/tau_er;
            }
            else {
                gI[i] = alpha(dt - dtI,tau_id,tau_ir,gI[i],hI[i]);
                dtI = dt;
                hI[i] = hI[i] - f/tau_ir;
            }
            is++;
            if (is == ns) {
                break;
            }
        }
    }
    if (extE) {
        if (noMid) {
            gEmid = alpha(tstep2-dtE,tau_ed,tau_er,gE[i],hE[i]);
            gE[i] = gEmid;
            dtE = tstep2;
        }
        gE[i] = alpha(tstep-dtE, tau_ed, tau_er, gE[i], hE[i]);
    }
    if (extI) {
        if (noMid) {
            gImid = alpha(tstep2-dtI,tau_id,tau_ir,gI[i],hI[i]);
            gI[i] = gImid;
            dtI = tstep2;
        }
        gI[i] = alpha(tstep-dtI, tau_id, tau_ir, gI[i], hI[i]);
    }
} 

void getCond4HH_old(std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &gEmid, std::vector<double> &gImid, std::vector<double> &tin, std::vector<size_b> &inID, double tau_er, double tau_ed, double tau_ir, double tau_id, std::vector<double> &fStrength, double tstep, size nt){
    const double conE = 1./(tau_er-tau_ed);
    const double conI = 1./(tau_ir-tau_id);
    double f,t,tstep2 = tstep/2;
    unsigned int i,j;
    for (j=0; j<tin.size(); j++) {
        f = fStrength[inID[j]];
        if (f > 0)
            for (i=ceil(tin[j]/tstep); i<nt; i++) {
                t = i*tstep;
                gE[i] = gE[i] + f*conE*(exp(-(t-tin[j])/tau_er)-exp(-(t-tin[j])/tau_ed));
                t = t + tstep2;
                gEmid[i] = gEmid[i] + f*conE*(exp(-(t-tin[j])/tau_er)-exp(-(t-tin[j])/tau_ed));
            }
        else
            for (i=ceil(tin[j]/tstep); i<nt; i++) {
                t = i*tstep;
                gI[i] = gI[i] - f*conI*(exp(-(t-tin[j])/tau_ir)-exp(-(t-tin[j])/tau_id));
                t = t + tstep2;
                gImid[i] = gImid[i] - f*conI*(exp(-(t-tin[j])/tau_ir)-exp(-(t-tin[j])/tau_id));
            }
    }
}

unsigned int RK4_HH_old(std::vector<double> &v, std::vector<double> &m, std::vector<double> &n,std::vector<double> &h, std::vector<double> &gE, std::vector<double> &gI, Neuron neuron, double pairs[], double tau_er, double tau_ed, double tau_ir, double tau_id, size nt, double tstep, double &tsp){
    unsigned int i, spikeCount = 0;
    double mk[4],nk[4],hk[4],fk[4];
    std::vector<double> gEmid(nt,0);
    std::vector<double> gImid(nt,0);
    double vT = pairs[8];
    double m_next, n_next, h_next, v_next, tstep2 = tstep/2;
    bool spiking = 0;
    getCond4HH_old(gE,gI,gEmid,gImid,neuron.tin,neuron.inID,tau_er,tau_ed,tau_ir,tau_id,neuron.fStrength,tstep,nt);
    for(i=1; i<nt; i++) {
        fk[0] = fk_HH(v[i-1],m[i-1],n[i-1],h[i-1],gE[i-1],gI[i-1],pairs);
        v_next = v[i-1] + fk[0]*tstep2;
        mk[0] = gatingm(v[i-1],m[i-1],vT);      m_next = m[i-1]+mk[0]*tstep2;
        nk[0] = gatingn(v[i-1],n[i-1],vT);      n_next = n[i-1]+nk[0]*tstep2; 
        hk[0] = gatingh(v[i-1],h[i-1],vT);      h_next = h[i-1]+hk[0]*tstep2;

        fk[1] = fk_HH(v_next,m_next,n_next,h_next, gEmid[i-1],gImid[i-1],pairs);
        v_next = v[i-1] + fk[1]*tstep2;
        mk[1] = gatingm(v_next,m_next,vT);      m_next = m[i-1]+mk[1]*tstep2;
        nk[1] = gatingn(v_next,n_next,vT);      n_next = n[i-1]+nk[1]*tstep2;
        hk[1] = gatingh(v_next,h_next,vT);      h_next = h[i-1]+hk[1]*tstep2;
        
        fk[2] = fk_HH(v_next,m_next,n_next,h_next, gEmid[i-1],gImid[i-1],pairs);
        v_next = v[i-1] + fk[2]*tstep2;
        mk[2] = gatingm(v_next,m_next,vT);      m_next = m[i-1]+mk[2]*tstep;
        nk[2] = gatingn(v_next,n_next,vT);      n_next = n[i-1]+nk[2]*tstep;
        hk[2] = gatingh(v_next,h_next,vT);      h_next = h[i-1]+hk[2]*tstep;
            
        fk[3] = fk_HH(v_next,m_next,n_next,h_next, gE[i],gI[i],pairs);
        mk[3] = gatingm(v_next,m_next,vT);   
        nk[3] = gatingn(v_next,n_next,vT);  
        hk[3] = gatingh(v_next,h_next,vT); 

        v[i] = v[i-1] + tstep/6*(fk[0] + 2*(fk[1]+fk[2]) + fk[3]);
        m[i] = m[i-1] + tstep/6*(mk[0] + 2*(mk[1]+mk[2]) + mk[3]);
        n[i] = n[i-1] + tstep/6*(nk[0] + 2*(nk[1]+nk[2]) + nk[3]);
        h[i] = h[i-1] + tstep/6*(hk[0] + 2*(hk[1]+hk[2]) + hk[3]);

        if (v[i]> neuron.vThres + (neuron.vThres-neuron.vReset)*2 && !spiking) {
            spikeCount = spikeCount + 1;
            if (v[i]<=v[i-1]) {
                tsp = (i-1)*tstep;
            }
            spiking = 1;
        } else if (v[i] < neuron.vThres) spiking = 0;
    }
    return spikeCount;
}

unsigned int RK4_IF(std::vector<double> &v, std::vector<double> &gE, std::vector<double> &gI, Neuron neuron, double *pairs, double tau_er, double tau_ed, double tau_ir, double tau_id, size nt, double tstep, bool cutoff, bool eif, double &tsp) {
    unsigned int i, spikeCount = 0;
    double fk[4], vThres, tstep2 = tstep/2;
    double gEr, gIr, gErmid, gIrmid, tstepsp, tstepsp2;
    std::vector<double> hE(nt,0), hI(nt,0);
    std::vector<double> gEmid(nt,0), gImid(nt,0);
    std::vector<std::vector<double>> tinBin(nt,std::vector<double>());
    std::vector<std::vector<double>> fBin(nt,std::vector<double>());
    if (eif) vThres = pairs[1];
    else vThres = pairs[9];
    getCond4IF(gE,gI,gEmid,gImid,hE,hI,tinBin,fBin,neuron.tin,neuron.inID,tau_er,tau_ed,tau_ir,tau_id,neuron.fStrength,tstep,nt);
    for(i=1; i<nt; i++) {
        fk[0] = fk_IF(               v[i-1],    gE[i-1],    gI[i-1],pairs,eif);
        fk[1] = fk_IF(v[i-1] + fk[0]*tstep2, gEmid[i-1], gImid[i-1],pairs,eif);
        fk[2] = fk_IF(v[i-1] + fk[1]*tstep2, gEmid[i-1], gImid[i-1],pairs,eif);
        fk[3] = fk_IF(v[i-1] + fk[2]*tstep2,      gE[i],      gI[i],pairs,eif);
        v[i] = v[i-1] + tstep/6*(fk[0] + 2*(fk[1]+fk[2]) + fk[3]);

        if (v[i]> vThres) {
            spikeCount = spikeCount + 1;
            //tsp = findTspRK(vThres,v[i-1],v[i],fk[0],fk_IF(v[i],gE[i],gI[i],pairs,eif));
            tsp = (vThres-v[i-1])/(v[i]-v[i-1])*tstep;
            tstepsp = tstep-tsp;
            tstepsp2 = tstepsp/2;
            getCond4t(gE[i-1],gI[i-1],hE[i-1],hI[i-1],gEr,gIr,gErmid,gIrmid,tinBin[i],fBin[i],tsp,(i-1)*tstep,tau_er,tau_ed,tau_ir,tau_id,tstep);

            fk[0] = fk_IF(               pairs[10],   gEr,   gIr,pairs,eif);
            fk[1] = fk_IF(pairs[10]+fk[0]*tstepsp2,gErmid,gIrmid,pairs,eif);
            fk[2] = fk_IF(pairs[10]+fk[1]*tstepsp2,gErmid,gIrmid,pairs,eif);
            fk[3] = fk_IF(pairs[10]+fk[2]*tstepsp2, gE[i], gI[i],pairs,eif);
            v[i] = pairs[10] + tstepsp/6*(fk[0] + 2*(fk[1] + fk[2]) + fk[3]);
            if (v[i] > vThres) std::cout << "use smaller tstep" << std::endl;
            if (cutoff) {
                tsp = (i-1) * tstep + tsp;
                break;
            }
        }
    }
    return spikeCount;
}	

unsigned int RK4_HH(std::vector<double> &v, std::vector<double> &m, std::vector<double> &n, std::vector<double> &h, std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &hE, std::vector<double> &hI, Neuron neuron, double pairs[], double tau_er, double tau_ed, double tau_ir, double tau_id, size nt, double tstep, double &tsp, bool insert, size vs, size &ve, size &is, double vStop){
    unsigned int i, j, spikeCount = 0;
    size vi, os, si = gE.size()-1;
    assert(gE.size() > 0);
    double mk[4],nk[4],hk[4],fk[4];
    double gEmid, gImid, t;
    double vT = pairs[8];
    double m_next, n_next, h_next, v_next, tstep2 = tstep/2;
    bool spiked,spiking = false;
    for(j=1; j<nt-vs; j++) {
        vi = vs + j;
        i = si + j;
        t = (vi-1)*tstep;
        getCond4HH(gE,gI,hE,hI,gEmid,gImid,neuron.tin,neuron.inID,is,neuron.tin.size(),tau_er,tau_ed,tau_ir,tau_id,neuron.fStrength,i,tstep,t,neuron.extE,neuron.extI);

        fk[0] = fk_HH(v[vi-1],m[i-1],n[i-1],h[i-1],gE[i-1],gI[i-1],pairs);
        v_next = v[vi-1] + fk[0]*tstep2;
        mk[0] = gatingm(v[vi-1],m[i-1],vT);      m_next = m[i-1]+mk[0]*tstep2;
        nk[0] = gatingn(v[vi-1],n[i-1],vT);      n_next = n[i-1]+nk[0]*tstep2;
        hk[0] = gatingh(v[vi-1],h[i-1],vT);      h_next = h[i-1]+hk[0]*tstep2;

        fk[1] = fk_HH(v_next,m_next,n_next,h_next, gEmid,gImid,pairs);
        v_next = v[vi-1] + fk[1]*tstep2;
        mk[1] = gatingm(v_next,m_next,vT);      m_next = m[i-1]+mk[1]*tstep2;
        nk[1] = gatingn(v_next,n_next,vT);      n_next = n[i-1]+nk[1]*tstep2;
        hk[1] = gatingh(v_next,h_next,vT);      h_next = h[i-1]+hk[1]*tstep2;
        
        fk[2] = fk_HH(v_next,m_next,n_next,h_next, gEmid,gImid,pairs);
        v_next = v[vi-1] + fk[2]*tstep;
        mk[2] = gatingm(v_next,m_next,vT);      m_next = m[i-1]+mk[2]*tstep;
        nk[2] = gatingn(v_next,n_next,vT);      n_next = n[i-1]+nk[2]*tstep;
        hk[2] = gatingh(v_next,h_next,vT);      h_next = h[i-1]+hk[2]*tstep;
            
        fk[3] = fk_HH(v_next,m_next,n_next,h_next, gE[i],gI[i],pairs);
        mk[3] = gatingm(v_next,m_next,vT);   
        nk[3] = gatingn(v_next,n_next,vT);  
        hk[3] = gatingh(v_next,h_next,vT); 

        v[vi] = v[vi-1] + tstep/6*(fk[0] + 2*(fk[1]+fk[2]) + fk[3]);
        m.push_back(m[i-1] + tstep/6*(mk[0] + 2*(mk[1]+mk[2]) + mk[3]));
        n.push_back(n[i-1] + tstep/6*(nk[0] + 2*(nk[1]+nk[2]) + nk[3]));
        h.push_back(h[i-1] + tstep/6*(hk[0] + 2*(hk[1]+hk[2]) + hk[3]));

        if (!spiking) {
            if (v[vi]> neuron.vThres) {
                spikeCount = spikeCount + 1;
                spiking = true;
                spiked = false;
            } 
        } else {
            if (v[vi]<=v[vi-1] && !spiked) {
                tsp = (vi-1)*tstep;
                std::cout << "spiked at " << tsp << std::endl;
                spiked = true;   
            }
            if (v[vi] < neuron.vThres) spiking = false;
        }
        if (insert) {
            if (spiked) {
                if (vi*tstep >= tsp + neuron.tref) {
                    ve = vi;
                    std::cout << "w/ spike " << v[vi] << std::endl;
                    return spikeCount;
                }
            } else {
                if (v[vi] < vStop) {
                    ve = vi; 
                    std::cout << "w/o spike " << v[vi] << std::endl;
                    return spikeCount;
                }
            }
        }
    }
    if (insert) {
        ve = nt-1;
        std::cout << "not back at the end" << std::endl;
    }
    return spikeCount;
}
#endif
