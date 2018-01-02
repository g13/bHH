#include <cmath>
#include <vector>
#include "typedefs.h"
#include "NeuronStruct.h"
#include "NeuronLibrary.h"
#include "singleRK4.h"
#include "singleRK2.h"

inline double interpTar(std::vector<double> &v, size i, double tstep, double t, double &tshift) {
    tshift = t - i*tstep;
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
inline void interpPSP(std::vector<double> &v, size vs, double ****PSP, double *vRange, size *idtRange,double *fRange, size nv, size ndt, size nf, double vTar, double dtTar, double fTar, double tshift, size ts, size tl, double *tmp) {
    size i[2],j[2],fi,fj,k,idt,jdt;
    double r[2],rf, base;
    getNear(vRange,nv,vTar,r[0],i[0],j[0]);
    getNear(idtRange,ndt,dtTar,r[1],i[1],j[1]);
    getNear(fRange,nf,fTar,rf,fi,fj);
    idt = idtRange[i[1]];
    jdt = idtRange[j[1]];
    
    for (k=0;k<tl;k++) {
        base = PSP[i[0]][i[1]][fi][idt+k];
        tmp[k] = base + r[0]*(PSP[j[0]][i[1]][fi][idt+k]-base) + r[1]*(PSP[i[0]][j[1]][fi][jdt+k]-base) + rf*(PSP[i[0]][i[1]][fj][idt+k]-base);
    }
    if (tshift!=1)
        for (k=1;k<tl;k++)
            v[vs+k] = v[vs+k] + tmp[k-1] + tshift*(tmp[k]-tmp[k-1]);
    else
        for (k=0;k<tl;k++)
            v[vs+k] = v[vs+k] + tmp[k];

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
    if (tshift!=1)
        for (k=1;k<tl;k++) {
            base = tmpK[k-1]*tmpPSP0[k-1]*PSP1[k-1];
            v[vs+k] = v[vs+k] + base + tshift*(tmpK[k]*tmpPSP0[k]*PSP1[k]-base);
        }
    else
        for (k=0;k<tl;k++)
            v[vs+k] = v[vs+k] + tmpK[k]*tmpPSP0[k]*PSP1[k];
}
inline void interpPSP0(std::vector<double> &v, size vs, double ***PSP0, double *vRange, double *fRange, size nv, size nf, double vTar, double fTar, double tshift, size ts, size tl, double *tmp) {
    size i[2],j[2],k;
    double r[2], base;
    getNear(vRange,nv,vTar,r[0],i[0],j[0]);
    getNear(fRange,nf,fTar,r[1],i[1],j[1]);
    
    for (k=0;k<tl;k++) {
        base = PSP0[i[0]][i[1]][ts+k];
        tmp[k] = base + r[0]*(PSP0[j[0]][i[1]][ts+k]-base) + r[1]*(PSP0[i[0]][j[1]][ts+k]-base);
    }
    if (tshift!=1)
        for (k=1;k<tl;k++)
            v[vs+k] = v[vs+k] + tmp[k-1] + tshift*(tmp[k]-tmp[k-1]);
    else
        for (k=0;k<tl;k++)
            v[vs+k] = v[vs+k] + tmp[k];
}
inline void interpVinit(std::vector<double> &v, size vs,  double **vLeak, double *vRange, size nv, double vTar, double tshift, size ts, size tl, double *tmp) {
    size i,j,k;
    double r;
    getNear(vRange,nv,vTar,r,i,j);
    if (tshift != 1) {
        for (k=0;k<tl;k++)
            tmp[k] = vLeak[i][ts+k] + r*(vLeak[j][ts+k]-vLeak[i][ts+k]);
        // [vs]---tmp[ts]-------[vs+1] ... [vs+l-1]---tmp[ts+l-1]
        //    <--->  tin   <-shift->                 <--->   ignored  
        for (k=1;k<tl;k++)
            v[vs+k] = tmp[k-1] + tshift*(tmp[k]-tmp[k-1]);
    } else for (k=0;k<tl;k++)
            v[vs+k] = vLeak[i][ts+k] + r*(vLeak[j][ts+k]-vLeak[i][ts+k]);
}
bool nearThreshold(Neuron neuron, NeuroLib neuroLib, std::vector<double> &v, std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &hE, std::vector<double> &hI, std::vector<double> &m, std::vector<double> &n, std::vector<double> &h, size vs, size nt, double tstep, double pairs[], double tau_er, double tau_ed, double tau_ir, double tau_id, double &tsp, double vinit, size &ve, size &ith, int rk, double vBack) {
    double conE = 1./(tau_er-tau_ed);
    double conI = 1./(tau_ir-tau_id);
    size i,j;
    unsigned int nc = 0;
    double f, lastCross = gE.size();
    double vT = pairs[9];
    gE.push_back(0);
    gI.push_back(0);
    hE.push_back(0);
    hI.push_back(0);
    double t = vs *tstep;
    for (j=ith+1; j>0; j--) {
        i = j-1;
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
    //if (rk==2)
    if (rk==4)
        ith++;
        nc = RK4_HH(v,m,n,h,gE,gI,hE,hI,neuron,pairs,tau_er,tau_ed,tau_ir,tau_id,nt,tstep,tsp,true,vs,ve,ith,vBack);
        ith--;
    if (nc) {
        assert(nc == 1);
        return true;
    } 
    else return false;
}
unsigned int bilinear_HH(std::vector<double> &v, NeuroLib neuroLib, Neuron neuron, double run_t, double ignore_t, double &tsp){
    std::string kType1, kType2, kType;
    double f1, f2, vTarget, dtTarget, tshift, shift;
    size ts, tl, vs, ve, i, j, k;
    size l0 = neuroLib.l0, nE = neuroLib.nE, nI = neuroLib.nI, nv = neuroLib.nv, nt0 = neuroLib.nt0-1, ndt = neuroLib.ndt;
    double *tmp = new double[nt0+1];
    double *tmpK = new double[nt0+1];
    double *tmpPSP = new double[nt0+1];
    double *currentPSP = new double[nt0+1];
    double tstep = neuroLib.tstep;
    size nt = static_cast<size>(run_t/tstep);
    double t0 = l0 - ignore_t;
    l0 = static_cast<size>(t0/tstep);
    unsigned int spikeCount = 0;

    
    double vinit = v[0];
    v.assign(nt+1,neuron.vReset);
    interpVinit(v,0,neuroLib.vLeak,neuroLib.vRange,nv,vinit,1,0,nt0,tmp);
    if (neuron.tin.size() > 0) {
        ve = static_cast<size>(neuron.tin[0]/tstep);
        for (i=0; i<neuron.tin.size(); i++) {
            // linear
            f2 = neuron.fStrength[neuron.inID[i]];
            vs = ve;
            vTarget = interpTar(v,vs,tstep,neuron.tin[i],tshift);
            shift = 1-tshift/tstep;

            if (vs+nt0 <= nt) 
                tl = nt0+1;
            else 
                tl = nt-vs+1;

            if (f2 > 0) {
                interpPSP0(v,vs, neuroLib.EPSP0, neuroLib.vRange, neuroLib.fE, nv, nE, vTarget, f2, shift,0,tl,currentPSP);
                kType2 = "1";
            } else  {
                interpPSP0(v,vs, neuroLib.IPSP0, neuroLib.vRange, neuroLib.fI, nv, nI, vTarget,-f2, shift,0,tl,currentPSP);
                kType2 = "0";
            }

            if (i<neuron.tin.size()-1) {
                ve = static_cast<size>(neuron.tin[i+1]/tstep);
            }
            else ve = nt+1;

            // bilinear 
            if (i!=0) {
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
                    if (dtTarget > t0) break;
                    dtTarget = dtTarget/tstep; // for interp along dtRange
                    if (run_t < neuron.tin[j] + t0)
                        tl = nt-vs+1;
                    else tl = static_cast<size>(ceil(l0 - dtTarget))+1;

                    f1 = neuron.fStrength[neuron.inID[j]];
                    if (f1 > 0) kType1 = "1";
                    else kType1 = "0";
                    kType = kType1 + kType2;
                    switch (std::stoi(kType,nullptr,10)) {
                        case 11: //EE
                            interpK(v, vs, neuroLib.kVEE, neuroLib.EPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fE, nv, ndt, nE, vTarget, dtTarget, f1, shift, tl, tmpK, tmpPSP);
                            break;
                        case 10: //EI
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
            }
            //check spike 
            for (k=vs;k<ve;k++) {
                if (k==0) continue;
                if (v[k] > neuron.vThres) {
                    tsp = (k-1+(neuron.vThres - v[k-1])/(v[k]-v[k-1])) * tstep;
                    neuron.tsp.push_back(tsp);
                    spikeCount = spikeCount + 1;
                    break;
                }
            }
        }
    }
    delete []tmp;
    delete []tmpK;
    delete []tmpPSP;
    delete []currentPSP;
    return spikeCount;
}

unsigned int linear_HH(std::vector<double> &v, std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &hE, std::vector<double> &hI, std::vector<double> &m, std::vector<double> &n, std::vector<double> &h, std::vector<size> &cross, NeuroLib neuroLib, Neuron neuron, double run_t, double ignore_t, double &tsp, double pairs[], double tau_er, double tau_ed, double tau_ir, double tau_id, double vCross, double vBack, double tref, int rk){
    double f, vTarget, dtTarget, tshift, shift;
    size ts, tl, vs, ve, vc, i, j, k, ith;
    size nE = neuroLib.nE, nI = neuroLib.nI, nv = neuroLib.nv, nt0 = neuroLib.nt0-1, ndt = neuroLib.ndt;
    double *tmp = new double[nt0+1];
    double *currentPSP = new double[nt0+1];
    double tstep = neuroLib.tstep;
    double tout;
    bool crossed;
    size ncross = 0;
    size nt = static_cast<size>(run_t/tstep);
    double t0 = tstep * nt0;
    double t1 = neuroLib.l0;
    size l1 = static_cast<size>(t1/tstep);
    unsigned int spikeCount = 0, spiked;

    
    double vinit = v[0];
    v.assign(nt+1,neuron.vReset);
    interpVinit(v,0,neuroLib.vLeak,neuroLib.vRange,nv,vinit,1,0,nt0,tmp);
    if (neuron.tin.size() > 0) {
        vs = static_cast<size>(neuron.tin[0]/tstep);
        ts = 0;
        for (i=0; i<neuron.tin.size(); i++) {
            // linear
            f = neuron.fStrength[neuron.inID[i]];
            vTarget = interpTar(v,vs,tstep,neuron.tin[i],tshift);
            shift = 1-tshift/tstep;

            if (vs+nt0 <= nt) 
                tl = nt0+1;
            else 
                tl = nt-vs+1;

            if (f > 0) {
                interpPSP0(v,vs, neuroLib.EPSP0, neuroLib.vRange, neuroLib.fE, nv, nE, vTarget, f, shift,0,tl,currentPSP);
            } else {
                interpPSP0(v,vs, neuroLib.IPSP0, neuroLib.vRange, neuroLib.fI, nv, nI, vTarget,-f, shift,0,tl,currentPSP);
            }

            
            if (i<neuron.tin.size()-1) {
                ve = static_cast<size>(neuron.tin[i+1]/tstep);
                if (ve > vs + nt0) {
                    vc = vs + nt0;
                } else {
                    vc = ve;
                }
            } else {
                vc = nt;
            }
            
            //deal consecutive nearThreshold until 
            crossed = false;
            ith = i;
            for (k=vs;k<vc;k++) {
                while (v[k] > vCross) {
                    crossed = true;
                    ncross = ncross +1;
                    // come in and out in multiples of tstep 
                    std::cout << "crossed at " << k*tstep << ", start ith " << ith << std::endl;
                    spiked = nearThreshold(neuron, neuroLib, v, gE, gI, hE, hI, m, n, h, k, nt+1, tstep, pairs, tau_er, tau_ed, tau_ir, tau_id, tsp, v[k], vs, ith, rk, vBack);
                    std::cout << "backed at " << vs*tstep << ", end ith " << ith << std::endl;
                    std::cout << "v[k] " << v[vs]  << " vCross " << vCross << " vBack " << vBack << std::endl;
                    cross.push_back(gE.size());
                    if (spiked){
                        spikeCount = spikeCount + 1;
                        neuron.tsp.push_back(tsp);
                    } 
                    if (vs + nt0 <= nt)
                        tl = nt0;
                    else tl = nt-vs+1;
                    interpVinit(v,vs,neuroLib.vLeak,neuroLib.vRange,nv,v[vs],1,0,tl,tmp);
                    vTarget = interpTar(v,vs,tstep,neuron.tin[j],tshift);
                    shift = 1-tshift/tstep;
                    if (vs < nt-1) {
                        for (i=ith+1;i>0; i--) {
                            j = i-1;
                            if (neuron.tin[j]/tstep + l1<= vs) {
                                break;
                            } else {
                                dtTarget = (vs - neuron.tin[j]/tstep);
                                f = neuron.fStrength[neuron.inID[j]];
                                if (run_t < neuron.tin[j] + t1)
                                    tl = nt-vs+1;
                                else tl = static_cast<size>(ceil(l1 - dtTarget))+1;
                                //std::cout << "t: " << neuron.tin[j] << " " << t1 << std::endl;
                                //std::cout << "it: " << dtTarget << " " << tl << " " << l1 << std::endl;
                                if (f > 0) {
                                    interpPSP(v,vs, neuroLib.EPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fE, nv, ndt, nE, vTarget, dtTarget, f, 1,0,tl,currentPSP);
                                } else {
                                    interpPSP(v,vs, neuroLib.IPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fI, nv, ndt, nI, vTarget, dtTarget,-f, 1,0,tl,currentPSP);
                                }
                            }
                        }
                        vc = vs + l1; 
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            if (ve < vc) {
                                vc = ve;
                            }
                        } else {
                            if (vc > nt)
                            vc = nt;
                        }

                        for (k=vs; k<vc; k++) {
                            if (v[k] > vCross) {
                                std::cout << " need to cross at " << k << std::endl;
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
