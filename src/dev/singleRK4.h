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
    // initialize g0, dtE,I (starting point of relaxation)
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
    // initialize bool for input(s) only exist in the second half of the time step;
    noMid = true; 

    if (is < ns) {
        // when there is input in the time step
        while (tin[is] < t + tstep) {
            dt = tin[is] - t;
            // if input is in the second half of time step
            // to calculate gMid, we only need to do one relaxation to gMid
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
                // gMid acquired.
                noMid = false;
            }
            // calculate relaxation toward input time, then jump h.
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
            // next input
            is++;
            if (is == ns) {
                break;
            }
        }
    }
    if (extE) {
        // if gMid is not acquired, meaning input time was in the first half of the time step or non-existed, calculate gMid.
        if (noMid) {
            gEmid = alpha(tstep2-dtE,tau_ed,tau_er,gE[i],hE[i]);
            gE[i] = gEmid;
            dtE = tstep2;
        }
        // final relaxation to next time step.
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

unsigned int RK4_HH(std::vector<double> &v, std::vector<double> &m, std::vector<double> &n, std::vector<double> &h, std::vector<double> &p, std::vector<double> &q, std::vector<double> &r, std::vector<double> &s, std::vector<double> &u, std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &hE, std::vector<double> &hI, Neuron neuron, double pairs[], bool type[], double tau_er, double tau_ed, double tau_ir, double tau_id, size nt, double tstep, std::vector<double> &tsp, bool insert, size vs, size &ve, size &is, double vStop){
    unsigned int i, j, spikeCount = 0;
    size vi, os, si = gE.size()-1;
    assert(gE.size() > 0);
    double mk[4],nk[4],hk[4],pk[4],qk[4],rk[4],uk[4],fk[4];
    double gEmid, gImid, t;
    double vT = pairs[8];
    double m_next = 0, n_next = 0, h_next = 0, v_next;
    double p_next = 0, q_next = 0, r_next = 0, s_next = 0, u_next = 0, tstep2 = tstep/2;
    bool spiked = false;
    bool spiking = false;
    double Vx = pairs[15];
    double tau_max = pairs[16];
    for(j=1; j<nt-vs; j++) {
        vi = vs + j;
        i = si + j;
        t = (vi-1)*tstep;
        getCond4HH(gE,gI,hE,hI,gEmid,gImid,neuron.tin,neuron.inID,is,neuron.tin.size(),tau_er,tau_ed,tau_ir,tau_id,neuron.fStrength,i,tstep,t,neuron.extE,neuron.extI);

        fk[0] = fk_HH(v[vi-1],m[i-1],n[i-1],h[i-1],p[i-1],q[i-1],r[i-1],s[i-1],u[i-1],gE[i-1],gI[i-1],pairs,type);
        if (type[0]) {
            mk[0] = gatingm(v[vi-1],m[i-1],vT);      
            hk[0] = gatingh(v[vi-1],h[i-1],vT);      
            m_next = m[i-1]+mk[0]*tstep2;
            h_next = h[i-1]+hk[0]*tstep2;
        }
        if (type[1]) {
            nk[0] = gatingn(v[vi-1],n[i-1],vT);      
            n_next = n[i-1]+nk[0]*tstep2;
        }
        if (type[2]) {
            pk[0] = gatingp(v[vi-1],p[i-1],tau_max);      
            p_next = p[i-1]+pk[0]*tstep2;
        }
        if (type[3]) {
            qk[0] = gatingq(v[vi-1],q[i-1]);
            rk[0] = gatingr(v[vi-1],r[i-1]);
            q_next = q[i-1]+qk[0]*tstep2;
            r_next = r[i-1]+rk[0]*tstep2;
        }
        if (type[4]) {
            uk[0] = gatingu(v[vi-1],u[i-1],Vx);
            u_next = u[i-1]+uk[0]*tstep2;
            v_next = v[vi-1] + fk[0]*tstep2;
            s_next = s_inf(v_next,Vx);
        } else {
            v_next = v[vi-1] + fk[0]*tstep2;
        }

        fk[1] = fk_HH(v_next,m_next,n_next,h_next,p_next,q_next,r_next,s_next,u_next,gEmid,gImid,pairs,type);
        if (type[0]) {
            mk[1] = gatingm(v_next,m_next,vT);
            hk[1] = gatingh(v_next,h_next,vT);
            m_next = m[i-1]+mk[1]*tstep2;
            h_next = h[i-1]+hk[1]*tstep2;
        }
        if (type[1]) {
            nk[1] = gatingn(v_next,n_next,vT);
            n_next = n[i-1]+nk[1]*tstep2;
        }
        if (type[2]) {
            pk[1] = gatingp(v_next,p_next,tau_max);      
            p_next = p[i-1]+pk[1]*tstep2;
        }
        if (type[3]) {
            qk[1] = gatingq(v_next,q_next);
            rk[1] = gatingr(v_next,r_next);
            q_next = q[i-1]+qk[1]*tstep2;
            r_next = r[i-1]+rk[1]*tstep2;
        }
        if (type[4]) {
            uk[1] = gatingu(v_next,u_next,Vx);
            u_next = u[i-1]+uk[1]*tstep2;
            v_next = v[vi-1] + fk[1]*tstep2;
            s_next = s_inf(v_next,Vx);
        } else {
            v_next = v[vi-1] + fk[1]*tstep2;
        }
        
        fk[2] = fk_HH(v_next,m_next,n_next,h_next,p_next,q_next,r_next,s_next,u_next,gEmid,gImid,pairs,type);
        if (type[0]) {
            mk[2] = gatingm(v_next,m_next,vT);
            hk[2] = gatingh(v_next,h_next,vT);
            m_next = m[i-1]+mk[2]*tstep;
            h_next = h[i-1]+hk[2]*tstep;
        }
        if (type[1]) {
            nk[2] = gatingn(v_next,n_next,vT);      
            n_next = n[i-1]+nk[2]*tstep;
        }
        if (type[2]) {
            pk[2] = gatingp(v_next,p_next,tau_max);      
            p_next = p[i-1]+pk[2]*tstep;
        }
        if (type[3]) {
            qk[2] = gatingq(v_next,q_next);
            rk[2] = gatingr(v_next,r_next);
            q_next = q[i-1]+qk[2]*tstep;
            r_next = r[i-1]+rk[2]*tstep;
        }
        if (type[4]) {
            uk[2] = gatingu(v_next,u_next,Vx);
            u_next = u[i-1]+uk[2]*tstep;
            v_next = v[vi-1]+fk[2]*tstep;
            s_next = s_inf(v_next,Vx);
        } else {
            v_next = v[vi-1]+fk[2]*tstep;
        }
            
        fk[3] = fk_HH(v_next,m_next,n_next,h_next,p_next,q_next,r_next,s_next,u_next,gE[i],gI[i],pairs,type);
        if (type[0]) {
            mk[3] = gatingm(v_next,m_next,vT);
            hk[3] = gatingh(v_next,h_next,vT); 
        }
        if (type[1]) {
            nk[3] = gatingn(v_next,n_next,vT);  
        }
        if (type[2]) {
            pk[3] = gatingp(v_next,p_next,tau_max);      
        }
        if (type[3]) {
            qk[3] = gatingq(v_next,q_next);
            rk[3] = gatingr(v_next,r_next);
        }
        if (type[4]) {
            uk[3] = gatingu(v_next,u_next,Vx);
        }

        v[vi] = v[vi-1] + tstep/6*(fk[0] + 2*(fk[1]+fk[2]) + fk[3]);
        if (type[0]) {
            m.push_back(m[i-1] + tstep/6*(mk[0] + 2*(mk[1]+mk[2]) + mk[3]));
            h.push_back(h[i-1] + tstep/6*(hk[0] + 2*(hk[1]+hk[2]) + hk[3]));
        } else { // pushback placeholder
            m.push_back(m_next);
            h.push_back(h_next);
        }
        if (type[1]){
            n.push_back(n[i-1] + tstep/6*(nk[0] + 2*(nk[1]+nk[2]) + nk[3]));
        } else {
            n.push_back(n_next);
        }
        if (type[2]) {
            p.push_back(p[i-1] + tstep/6*(pk[0] + 2*(pk[1]+pk[2]) + pk[3]));
        } else {
            p.push_back(p_next);
        }
        if (type[3]) {
            q.push_back(q[i-1] + tstep/6*(qk[0] + 2*(qk[1]+qk[2]) + qk[3]));
            r.push_back(r[i-1] + tstep/6*(rk[0] + 2*(rk[1]+rk[2]) + rk[3]));
        } else {
            q.push_back(q_next);
            r.push_back(r_next);
        }
        if (type[4]) {
            u.push_back(u[i-1] + tstep/6*(uk[0] + 2*(uk[1]+uk[2]) + uk[3]));
            s.push_back(s_inf(v[vi],Vx));
        } else {
            u.push_back(u_next);
            s.push_back(s_next);
        }
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
