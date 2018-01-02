#include <cmath>
#include <vector>
#include "typedefs.h"
#include "NeuronStruct.h"
#include "NeuronLibrary.h"
#include "singleRK4.h"
#include "singleRK2.h"
namespace sB{
    const bool debug = true;
    const bool debug2 = false;
}
inline double interpTar(std::vector<double> &v, size i, double tstep, double t, double &tshift) {
    tshift = t/tstep - i;
    if (tshift==0)
        return v[i];
    else
        return v[i] + tshift*(v[i+1]-v[i]);
}
template<typename T>
inline void getNear(T *range, size n, double target, double &ratio, size &istart, size &jnext) {
    size i;
    if (target < range[0]) {
        istart = 1; jnext = 0;
        ratio = (target-range[1])/(range[0]-range[1]);
    } else if (target >= range[n-1]) {
        istart = n-2; jnext = n-1;
        ratio = (target-range[istart])/(range[jnext]-range[istart]);
    } else {
        for (i=0;i<n-1;i++) {
            if (target >= range[i] && target < range[i+1] ) {
                istart = i; jnext = i+1;
                ratio = (target-range[istart])/(range[jnext]-range[istart]);
                break;
            }
        }
    }
}
inline void interpPSP(std::vector<double> &v, size vs, double ****PSP, double *vRange, size *idtRange,double *fRange, size nv, size ndt, size nf, double vTar, double dtTar, double fTar, double *tmp, size limit) {
    size i[2],j[2],fi,fj,k,idt,jdt;
    double r[2],rf, base;
    getNear(vRange,nv,vTar,r[0],i[0],j[0]);
    getNear(idtRange,ndt,dtTar,r[1],i[1],j[1]);
    getNear(fRange,nf,fTar,rf,fi,fj);
    idt = idtRange[i[1]];
    jdt = idtRange[j[1]];
    
    size tl = limit - jdt + 1;

    for (k=1;k<tl;k++) {
        base = PSP[i[0]][i[1]][fi][idt+k];
        tmp[k] = base + r[0]*(PSP[j[0]][i[1]][fi][idt+k]-base) + r[1]*(PSP[i[0]][j[1]][fi][jdt+k]-base) + rf*(PSP[i[0]][i[1]][fj][idt+k]-base);
    }
    for (k=1;k<tl;k++) {
        v[vs+k] += tmp[k];
    }
}

inline void interpK(std::vector<double> &v, size vs, double ***kV, double ****PSP0, double *PSP1, double *vRange, size *idtRange, double *fRange, size nv, size ndt, size nf,  double vTar, double dtTar, double fTar, double tshift, size tl, double *tmpK, double *tmpPSP0) {
    size i[2],j[2],fi,fj,k,idt,jdt;
    double r[2],rf, base;
    getNear(vRange,nv,vTar,r[0],i[0],j[0]);
    getNear(idtRange,ndt,dtTar,r[1],i[1],j[1]);
    getNear(fRange,nf,fTar,rf,fi,fj);
    idt = idtRange[i[1]];
    jdt = idtRange[j[1]];
    
    for (k=0;k<tl;k++) {
        base = kV[i[0]][i[1]][idt+k];
        tmpK[k] = base + r[0]*(kV[j[0]][i[1]][idt+k]-base) + r[1]*(kV[i[0]][j[1]][jdt+k]-base);
        base = PSP0[i[0]][i[1]][fi][idt+k];
        tmpPSP0[k] = base + r[0]*(PSP0[j[0]][i[1]][fi][idt+k]-base) + r[1]*(PSP0[i[0]][j[1]][fi][jdt+k]-base) + rf*(PSP0[i[0]][i[1]][fj][idt+k]-base);
    }
    if (tshift!=1) {
        for (k=1;k<tl;k++) {
            base = tmpK[k-1]*tmpPSP0[k-1]*PSP1[k-1];
            v[vs+k] = v[vs+k] + base + tshift*(tmpK[k]*tmpPSP0[k]*PSP1[k]-base);
        }
    } else {
        for (k=1;k<tl;k++) {
            v[vs+k] = v[vs+k] + tmpK[k]*tmpPSP0[k]*PSP1[k];
        }
    }
}

inline void interpPSP0(std::vector<double> &v, size vs, double ***PSP0, double *vRange, double *fRange, size nv, size nf, double vTar, double fTar, double tshift, size tl, double *tmp) {
    size i[2],j[2],k;
    double r[2], base;
    getNear(vRange,nv,vTar,r[0],i[0],j[0]);
    getNear(fRange,nf,fTar,r[1],i[1],j[1]);
    
    //std::cout << "  " << r[0] << ", " << i[0] << ", " << j[0] << std::endl;
    for (k=1;k<tl;k++) {
        base = PSP0[i[0]][i[1]][k];
        tmp[k] = base + r[0]*(PSP0[j[0]][i[1]][k]-base) + r[1]*(PSP0[i[0]][j[1]][k]-base);
    }
    if (tshift==1) {
        for (k=1;k<tl;k++) {
            v[vs+k] += tmp[k];
        }
    } else {
        for (k=1;k<tl;k++) {
            v[vs+k] += tmp[k-1] + tshift*(tmp[k]-tmp[k-1]);
        }
    }
}
inline void interpVinit(std::vector<double> &v, size vs,  double **vLeak, double *vRange, size nv, double vTar, size tl, double *tmp) {
    size i,j,k;
    double r;
    getNear(vRange,nv,vTar,r,i,j);
    for (k=0;k<tl;k++)
        v[vs+k] = vLeak[i][k] + r*(vLeak[j][k]-vLeak[i][k]);
}
bool nearThreshold(Neuron &neuron, NeuroLib &neuroLib, std::vector<double> &v, std::vector<double> &crossv, std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &hE, std::vector<double> &hI, std::vector<double> &m, std::vector<double> &n, std::vector<double> &h, size vs, size nt, double tstep, double pairs[], double tau_er, double tau_ed, double tau_ir, double tau_id, std::vector<double> &tsp, double vinit, size &ve, size &ith, int rk, double vBack) {
    double conE = 1./(tau_er-tau_ed);
    double conI = 1./(tau_ir-tau_id);
    size i,j;
    unsigned int nc = 0;
    double f;
    size endCross, lCross, lastCross = gE.size();
    double vT = pairs[8];
    gE.push_back(0);
    gI.push_back(0);
    hE.push_back(0);
    hI.push_back(0);
    double t = vs *tstep;
    for (j=ith+1; j>0; j--) {
        i = j-1;
        //std::cout << "eeeee " << i << std::endl;
        if (neuron.tin[i] + neuroLib.nt0*tstep > t) {
            f = neuron.fStrength[neuron.inID[i]];
            if (f > 0) {
                gE[lastCross] = gE[lastCross] + f*conE*(exp(-(t-neuron.tin[i])/tau_er)-exp(-(t-neuron.tin[i])/tau_ed));
                hE[lastCross] = hE[lastCross] + f/tau_er*exp(-(t-neuron.tin[i])/tau_er);
            } else {
                gI[lastCross] = gI[lastCross] - f*conI*(exp(-(t-neuron.tin[i])/tau_ir)-exp(-(t-neuron.tin[i])/tau_id));
                hI[lastCross] = hI[lastCross] - f/tau_ir*exp(-(t-neuron.tin[i])/tau_ir);
            }
        } else break;
    }
    m.push_back(m_inf(vinit,vT));
    n.push_back(n_inf(vinit,vT));
    h.push_back(h_inf(vinit,vT));
    crossv.push_back(v[vs]);
    //if (rk==2)
    ith++;
    if (rk==4) {
        nc = RK4_HH(crossv,m,n,h,gE,gI,hE,hI,neuron,pairs,tau_er,tau_ed,tau_ir,tau_id,nt,tstep,tsp,true,lastCross,endCross,ith,vBack,t);
    }
    lCross = endCross - lastCross;
    for (j=1;j<=lCross;j++) {
        v[vs+j] = crossv[lastCross+j];
    }
    ve = vs + lCross;
    if (ve > nt) {
        std::cout << "ve " << ve << std::endl;
        assert(ve < nt);
    }
    ith--;
    if (nc) {
        assert(nc == 1);
        return true;
    } 
    else return false;
}

unsigned int bilinear_HH(std::vector<double> &v, std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &hE, std::vector<double> &hI, std::vector<double> &m, std::vector<double> &n, std::vector<double> &h, std::vector<double> &crossv, NeuroLib &neuroLib, Neuron &neuron, double run_t, double ignore_t, std::vector<double> &tsp, double pairs[], double tau_er, double tau_ed, double tau_ir, double tau_id, double vCross, double vBack, double tref, int rk){
    double f1, f2, vTarget, dtTarget, tshift, shift;
    size tl, vs, ve, vc, i, ii, j, k, ith;
    double tstep = neuroLib.tstep;
    size nE = neuroLib.nE;
    size nI = neuroLib.nI;
    size nv = neuroLib.nv;
    size ndt = neuroLib.ndt;
    size nt0 = neuroLib.nt0;
    size l0 = nt0 - 1;
    size t0 = l0 * tstep;

    size nt = neuroLib.nt;
    size l1 = nt - 1;
    size t1 = l1 * tstep;

    size tb = neuroLib.l0-ignore_t;
    size lb = static_cast<size>(tb/tstep);
    size nb = lb + 1;

    size run_lt = static_cast<size>(run_t/tstep);
    size run_nt = run_lt + 1;
    unsigned int kType;
    double *tmp = new double[nt0];
    double *tmpK = new double[nt0];
    double *tmpPSP = new double[nt0];
    double *currentPSP = new double[nt0];
    currentPSP[0] = 0;
    bool crossed; 
    size ncross = 0;
    
    size idt = ndt-1;
    size lastInterval = nt - neuroLib.idtRange[idt];
    while (lastInterval < 0) {
        idt--;
        lastInterval = nt - neuroLib.idtRange[idt];
    }
    double lastIntervalt = (lastInterval-1)*tstep;
    
    size limit;
    unsigned int spikeCount = 0, spiked;
    
    double vinit = v[0];
    v.assign(run_nt,neuron.vReset);
    if (l0 > run_lt) {
        tl = run_nt;
    } else {
        tl = nt0;
    }
    interpVinit(v,0,neuroLib.vLeak,neuroLib.vRange,nv,vinit,tl,tmp);
    if (neuron.tin.size() > 0) {
        vs = static_cast<size>(neuron.tin[0]/tstep);
        for (i=0; i<neuron.tin.size(); i++) {
            // linear
            f2 = neuron.fStrength[neuron.inID[i]];
            vTarget = interpTar(v,vs,tstep,neuron.tin[i],tshift);
            shift = 1-tshift;

            if (vs+l0 > run_lt) {
                tl = run_nt-vs;
            } else {
                tl = nt0;
            }

            if (f2 > 0) {
                interpPSP0(v,vs, neuroLib.EPSP0, neuroLib.vRange, neuroLib.fE, nv, nE, vTarget, f2, shift,tl,currentPSP);
                kType = 1;
            } else  {
                interpPSP0(v,vs, neuroLib.IPSP0, neuroLib.vRange, neuroLib.fI, nv, nI, vTarget,-f2, shift,tl,currentPSP);
                kType = 0;
            }

            if (i<neuron.tin.size()-1) {
                ve = static_cast<size>(neuron.tin[i+1]/tstep);
                if (ve + 1 > vs + nt0) {
                    vc = vs + nt0;
                } else {
                    vc = ve + 1;
                }
            } else {
                vc = run_nt;
            }
            // bilinear 
                //      <-------------------------l0--------------->
                //  tin[j]-------------------tin[i]------------tin[j]+l0------tin[i]+dur----tin[i]+dur0
                //     *------------------------*------------------*---------------*-------------*----
                // tstep --|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
                //     <---> 1-shift[j]      <--> tshift[i]
                // k0dt|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
                // kpdt ----------NA------------|-----|-----|-----|-----|-----|-----|-----|-----|-----|
                //  kJI < ----------0----------> <---------------1----------------> <------NA---->
                //  biV <-----------0----------> <--------1-------> <------------NA-------------->
                //      <------dtTarget--------> <-------tl-------> 
                //                                                  <-------------->
                //                                                 interp space for dt 
            for (k=i; k>0; k--) {
                j = k-1; // prevent negative for unsigned int iteration number. 
                dtTarget = neuron.tin[i]-neuron.tin[j];
                if (dtTarget > tb) break;
                dtTarget = dtTarget/tstep; // for interp along idtRange
                if (neuron.tin[j] + tb > run_t)
                    tl = run_nt - vs;
                else tl = static_cast<size>(floor(nb - dtTarget - shift));

                f1 = neuron.fStrength[neuron.inID[j]];
                if (f1 > 0) kType = kType + 2;
                switch (kType) {
                    case 3: //EE
                        interpK(v, vs, neuroLib.kVEE, neuroLib.EPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fE, nv, ndt, nE, vTarget, dtTarget, f1, shift, tl, tmpK, tmpPSP);
                        break;
                    case 2: //EI
                        interpK(v, vs, neuroLib.kVEI, neuroLib.EPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fE, nv, ndt, nI, vTarget, dtTarget, f1, shift, tl, tmpK, tmpPSP);
                        break;
                    case 1:  //IE
                        interpK(v, vs, neuroLib.kVIE, neuroLib.IPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fI, nv, ndt, nE, vTarget, dtTarget, -f1, shift, tl, tmpK, tmpPSP);
                        break;
                    case 0:  //II
                        interpK(v, vs, neuroLib.kVII, neuroLib.IPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fI, nv, ndt, nI, vTarget, dtTarget, -f1, shift, tl, tmpK, tmpPSP);
                        break;
                }
            } 
            //check spike 
            crossed = false;
            ith = i;
            for (k=vs;k<vc;k++) {
                while(v[k] > vCross) {
                    std::cout << "k " << k << std::endl;
                    assert(k<vc);
                    crossed = true;
                    ncross = ncross +1;
                    if (sB::debug) {
                        std::cout << "crossed at " << k*tstep << ", start ith " << ith << std::endl;
                    }
                    spiked = nearThreshold(neuron, neuroLib, v, crossv, gE, gI, hE, hI, m, n, h, k, run_nt, tstep, pairs, tau_er, tau_ed, tau_ir, tau_id, tsp, v[k], vs, ith, rk, vBack);
                    if (sB::debug) {
                        std::cout << "backed at " << vs*tstep << ", end ith " << ith << std::endl;
                        std::cout << "v[k] " << v[vs]  << " vCross " << vCross << " vBack " << vBack << std::endl;
                    }
                    if (spiked){
                        spikeCount = spikeCount + 1;
                        neuron.tsp.push_back(tsp.back());
                    } 
                    if (vs + l0 <= run_lt)
                        tl = nt0;
                    else tl = run_nt-vs;
                    vTarget = v[vs];
                    interpVinit(v,vs,neuroLib.vLeak,neuroLib.vRange,nv,vTarget,tl,tmp);
                    if (vs < run_lt) {
                        for (ii=ith+1; ii>0; ii--) {
                            i = ii-1;
                            dtTarget = vs - neuron.tin[i]/tstep;
                            if (l1 - lastInterval <= dtTarget) {
                                break;
                            }
                            f2 = neuron.fStrength[neuron.inID[i]];

                            if (neuron.tin[i] + t1 > run_t) {
                                limit = run_lt - vs;
                            }
                            else {
                                limit = l1;
                            }

                            if (f2 > 0) {
                                interpPSP(v,vs, neuroLib.EPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fE, nv, ndt, nE, vTarget, dtTarget, f2,currentPSP, limit);
                                kType = 1;
                            } else {
                                interpPSP(v,vs, neuroLib.IPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fI, nv, ndt, nI, vTarget, dtTarget,-f2,currentPSP, limit);
                                kType = 0;
                            }
                            for (k=i; k>0; k--) {
                                j = k-1; // prevent negative for unsigned int iteration number. 
                                dtTarget = vs - neuron.tin[j]/tstep;
                                if (dtTarget > lb) {
                                    break;
                                }
                                if (neuron.tin[j] + tb > run_t)
                                    tl = run_nt - vs;
                                else tl = static_cast<size>(floor(nb - dtTarget));

                                f1 = neuron.fStrength[neuron.inID[j]];
                                if (f1 > 0) kType = kType + 2;
                                switch (kType) {
                                    case 3: //EE
                                        interpK(v, vs, neuroLib.kVEE, neuroLib.EPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fE, nv, ndt, nE, vTarget, dtTarget, f1, 1, tl, tmpK, tmpPSP);
                                        break;
                                    case 2: //EI
                                        interpK(v, vs, neuroLib.kVEI, neuroLib.EPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fE, nv, ndt, nI, vTarget, dtTarget, f1, 1, tl, tmpK, tmpPSP);
                                        break;
                                    case 1:  //IE
                                        interpK(v, vs, neuroLib.kVIE, neuroLib.IPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fI, nv, ndt, nE, vTarget, dtTarget, -f1, 1, tl, tmpK, tmpPSP);
                                        break;
                                    case 0:  //II
                                        interpK(v, vs, neuroLib.kVII, neuroLib.IPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fI, nv, ndt, nI, vTarget, dtTarget, -f1, 1, tl, tmpK, tmpPSP);
                                        break;
                                }
                            } 
                        }
                        vc = vs + l1; 
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            if (ve + 1 < vc) {
                                vc = ve + 1;
                            }
                        } else {
                            if (vc > run_lt) {
                                vc = run_nt;
                            }
                        }
                        k = vc-1;
                        for (j=vs; j<vc; j++) {
                            if (v[j] > vCross) {
                                k = j;
                                if (sB::debug) {
                                    std::cout << "cross again at " << k << std::endl;
                                }
                                break;
                            }
                        }
                    } else {
                        break;
                    }
                }
                if (crossed) {
                    i = ith;
                    break;
                }
            }
            vs = ve;
        }
    }
    delete []tmp;
    delete []tmpK;
    delete []tmpPSP;
    delete []currentPSP;
    std::cout << "crossed " << ncross << " times" << std::endl;
    return spikeCount;
}

unsigned int linear_HH(std::vector<double> &v, std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &hE, std::vector<double> &hI, std::vector<double> &m, std::vector<double> &n, std::vector<double> &h, std::vector<double> &crossv, NeuroLib &neuroLib, Neuron &neuron, double run_t, double ignore_t, std::vector<double> &tsp, double pairs[], double tau_er, double tau_ed, double tau_ir, double tau_id, double vCross, double vBack, double tref, int rk){
    double f, vTarget, dtTarget, tshift, shift;
    size tl, vs, ve, vc, i, j, k, ith;
    double tstep = neuroLib.tstep;
    size nE = neuroLib.nE;
    size nI = neuroLib.nI;
    size nv = neuroLib.nv;
    size nt0 = neuroLib.nt0;
    size ndt = neuroLib.ndt;
    size l0 = nt0 - 1;
    size t0 = l0 * tstep;

    size nt = neuroLib.nt;
    size l1 = nt - 1;
    size t1 = l1 * tstep;

    size run_lt = static_cast<size>(run_t/tstep);
    size run_nt = run_lt + 1;

    double *tmp = new double[nt0+1];
    double *currentPSP = new double[nt0+1];
    currentPSP[0] = 0;
    bool crossed;
    size ncross = 0;

    size idt = ndt-1;
    size lastInterval = nt - neuroLib.idtRange[idt];
    while (lastInterval < 0) {
        idt--;
        lastInterval = nt - neuroLib.idtRange[idt];
    }
    double lastIntervalt = (lastInterval-1)*tstep;

    size limit;
    unsigned int spikeCount = 0, spiked;
    
    double vinit = v[0];
    v.assign(run_nt,neuron.vReset);
    if (l0 > run_lt) {
        tl = run_nt;
    } else {
        tl = nt0;
    }
    interpVinit(v,0,neuroLib.vLeak,neuroLib.vRange,nv,vinit,tl,tmp);
    if (neuron.tin.size() > 0) {
        vs = static_cast<size>(neuron.tin[0]/tstep);
        for (i=0; i<neuron.tin.size(); i++) {
            // linear
            f = neuron.fStrength[neuron.inID[i]];
            vTarget = interpTar(v,vs,tstep,neuron.tin[i],tshift);
            shift = 1-tshift;

            if (vs+l0 > run_lt) {
                tl = run_nt-vs;
            } else {
                tl = nt0;
            }

            if (f > 0) {
                interpPSP0(v,vs, neuroLib.EPSP0, neuroLib.vRange, neuroLib.fE, nv, nE, vTarget, f, shift,tl,currentPSP);
            } else {
                interpPSP0(v,vs, neuroLib.IPSP0, neuroLib.vRange, neuroLib.fI, nv, nI, vTarget,-f, shift,tl,currentPSP);
            }

            
            if (i<neuron.tin.size()-1) {
                ve = static_cast<size>(neuron.tin[i+1]/tstep);
                if (ve + 1 > vs + nt0) {
                    vc = vs + nt0;
                } else {
                    vc = ve + 1;
                }
            } else {
                vc = run_nt;
            }
            
            //deal consecutive nearThreshold until 
            crossed = false;
            ith = i;
            for (k=vs;k<vc;k++) {
                while (v[k] > vCross) {
                    std::cout << "k " << k << std::endl;
                    crossed = true;
                    ncross = ncross +1;
                    // come in and out in multiples of tstep 
                    if (sB::debug) {
                        std::cout << "crossed at " << k*tstep << ", start ith " << ith << std::endl;
                    }
                    spiked = nearThreshold(neuron, neuroLib, v, crossv, gE, gI, hE, hI, m, n, h, k, run_nt, tstep, pairs, tau_er, tau_ed, tau_ir, tau_id, tsp, v[k], vs, ith, rk, vBack);
                    if (sB::debug) {
                        std::cout << "backed at " << vs*tstep << ", end ith " << ith << std::endl;
                        std::cout << "v[k] " << v[vs]  << " vCross " << vCross << " vBack " << vBack << std::endl;
                    }
                    //cross.push_back(gE.size());
                    if (spiked){
                        spikeCount = spikeCount + 1;
                        neuron.tsp.push_back(tsp.back());
                    } 
                    if (vs + l0 <= run_lt)
                        tl = nt0;
                    else tl = run_nt-vs;
                    vTarget = v[vs];
                    interpVinit(v,vs,neuroLib.vLeak,neuroLib.vRange,nv,vTarget,tl,tmp);
                    if (vs < run_lt) {
                        for (i=ith+1;i>0; i--) {
                            j = i-1;
                            if (neuron.tin[j]/tstep + l1 - lastInterval <= vs) {
                                break;
                            }
                            dtTarget = (vs - neuron.tin[j]/tstep);
                            f = neuron.fStrength[neuron.inID[j]];

                            if (neuron.tin[j] + t1 > run_t) {
                                limit = run_lt - vs;
                            } else {
                                limit = l1;
                            }
                            if (f > 0) {
                                interpPSP(v,vs, neuroLib.EPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fE, nv, ndt, nE, vTarget, dtTarget, f, currentPSP,limit);
                            } else {
                                interpPSP(v,vs, neuroLib.IPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fI, nv, ndt, nI, vTarget, dtTarget,-f,currentPSP,limit);
                            }
                        }
                        vc = vs + l1; 
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            if (ve < vc) {
                                vc = ve;
                            }
                        } else {
                            vc = run_nt;
                        }
                        k = vc-1;
                        for (j=vs; j<vc; j++) {
                            if (v[j] > vCross) {
                                k = j;
                                if (sB::debug) {
                                    std::cout << "cross again at " << k << std::endl;
                                }
                                break;
                            }
                        }
                    } else break;
                }
                if (crossed) {
                    i = ith;
                    break;
                }
            }
            vs = ve;
        }
    }
    delete []tmp;
    delete []currentPSP;
    std::cout << "crossed " << ncross << " times" << std::endl;
    return spikeCount;
}
