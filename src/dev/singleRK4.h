#ifndef SINGLE_RK4_H
#define SINGLE_RK4_H

#include "NeuronStruct.h"
#include "math.h"
#include <vector>
#include <iostream>
#include "channels.h"
#include "typedefs.h"

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

unsigned int RK4_HH(std::vector<double> &v, std::vector<double> &m, std::vector<double> &n, std::vector<double> &h, std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &hE, std::vector<double> &hI, Neuron neuron, double pairs[], double tau_er, double tau_ed, double tau_ir, double tau_id, size nt, double tstep, std::vector<double> &tsp, bool insert, size vs, size &ve, size &is, double vStop){
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
        if (insert) {
            //std::cout << "fk: " << fk[0] << ", " << fk[1] << ", " << fk[2] << ", " << fk[3] << std::endl;
        }

        if (!spiking) {
            if (v[vi]> neuron.vThres) {
                spikeCount = spikeCount + 1;
                spiking = true;
                spiked = false;
            } 
        } else {
            if (v[vi]<=v[vi-1] && !spiked) {
                tsp.push_back((vi-1)*tstep);
                std::cout << "spiked at " << tsp.back() << std::endl;
                spiked = true;   
            }
            if (v[vi] < neuron.vThres) spiking = false;
        }
        if (insert) {
            if (spiked) {
                if ((vi*tstep >= tsp.back() + neuron.tref) && v[vi] < vStop) {
                    ve = vi;
                    std::cout << "w/ spike " << v[vi] << ", passed " << j*tstep << "ms" << std::endl;
                    return spikeCount;
                }
            } else {
                if (v[vi] < vStop) {
                    ve = vi; 
                    std::cout << "w/o spike " << v[vi] << ", passed " << j*tstep << "ms" << std::endl;
                    return spikeCount;
                }
            }
        }
    }
    if (insert) {
        ve = nt-1;
        std::cout << "not back at the end, passed " << j*tstep << "ms" << std::endl;
    }
    return spikeCount;
}

#endif
