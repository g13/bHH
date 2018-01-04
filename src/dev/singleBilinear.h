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

inline void interpkV(std::vector<double> &v, size vs, double ******kV, int iSyn, int jSyn, double *vRange, size *idtRange, size nv, size ndt,  double vTar, double dtTar0, double dtTar1, double tshift, size tl, double *tmpkV) {
    size idt,jdt;
    size fi,fj;
    size idt0,jdt0;
    size idt1,jdt1;
    size iv, jv;
    double rv, rdt, rf, base;
    getNear(vRange,nv,vTar,rv,iv,jv);
    getNear(idtRange,ndt,dtTar0,rdt,idt,jdt);
    double dtTar = dtTar1 - dtTar0;
    if (dtTar == 0) {
        size iidt = idtRange[idt];
        size jjdt = jdtRange[jdt];
        for (size i=0;i<tl;i++) {
            base = kV[iv][idt][iSyn][jSyn][0][iidt+i];
            tmpK[i] = base + rv*(kV[jv][idt][iSyn][jSyn][0][iidt+i]-base) + rdt*(kV[iv][jdt][iSyn][jSyn][0][jjdt+i]-base);
        }
    } else {
        getNear(idtRange,ndt,idtRange[idt] + dtTar,rdt0,idt0,jdt0);
        getNear(idtRange,ndt,idtRange[jdt] + dtTar,rdt1,idt1,jdt1);
        size iidt0 = idtRange[idt0];
        size jjdt0 = idtRange[jdt0];
        size iidt1 = idtRange[idt1];
        size jjdt1 = idtRange[jdt1];
        for (size i=0;i<tl;i++) {
            double base0 = kV[iv][idt][iSyn][jSyn][idt0][iidt0+i];
            base = base0 + rdt0*(kV[iv][idt][iSyn][jSyn][jdt0][jjdt0+i] - base0);
            
            double dv0 = kV[jv][idt][iSyn][jSyn][idt0][iidt0+i];
            dv = dv0 + rdt0*(kV[jv][idt][iSyn][jSyn][jdt0][jjdt0+i] - dv0);

            double dt0 = kV[iv][jdt][iSyn][jSyn][idt1][iidt1+i];
            dt = dt0 + rdt1*(kV[iv][jdt][iSyn][jSyn][jdt1][jjdt1+i] - dt0);

            tmpkV[i] = base + rv*(dv-base) + rdt*(dt-base);
        }
    }

    if (tshift!=1) {
        for (size i=1;i<tl;i++) {
            v[vs+i] += tmpkV[i-1] + tshift*(tmpkV[i]-tmpkV[i-1]);
        }
    } else {
        for (size i=0;i<tl;i++) {
            v[vs+i] += tmpkV[i];
        }
    }
}

inline void interpK(std::vector<double> &v, size vs, double ****k, double ****PSP0, double *PSP1, double *vRange, size *idtRange, double *fRange, size nv, size ndt, size nf,  double vTar, double dtTar0, double dtTar1, double fTar, double tshift, size tl, double *tmpK, double *tmpPSP0) {
    size idt,jdt;
    size fi,fj;
    size idt0,jdt0;
    size idt1,jdt1;
    size iv, jv;
    double rv, rdt, rf, base;
    getNear(vRange,nv,vTar,rv,iv,jv);
    getNear(idtRange,ndt,dtTar0,rdt,idt,jdt);
    getNear(fRange,nf,fTar,rf,fi,fj);
    double dtTar = dtTar1 - dtTar0;
    if (dtTar == 0) {
        size iidt = idtRange[idt];
        size jjdt = jdtRange[jdt];
        for (size i=0;i<tl;i++) {
            base = k[iv][idt][iidt+i];
            tmpK[i] = base + rv*(kV[jv][idt][0][iidt+i]-base) + rdt*(kV[iv][jdt][0][jjdt+i]-base);
        }
    } else {
        getNear(idtRange,ndt,idtRange[idt] + dtTar,rdt0,idt0,jdt0);
        getNear(idtRange,ndt,idtRange[jdt] + dtTar,rdt1,idt1,jdt1);
        size iidt0 = idtRange[idt0];
        size jjdt0 = idtRange[jdt0];
        size iidt1 = idtRange[idt1];
        size jjdt1 = idtRange[jdt1];
        for (size i=0;i<tl;i++) {
            double base0 = k[iv][idt][idt0][iidt0+i];
            base = base0 + rdt0*(k[iv][idt][jdt0][jjdt0+i] - base0);
            
            double dv0 = k[jv][idt][idt0][iidt0+i];
            dv = dv0 + rdt0*(k[jv][idt][jdt0][jjdt0+i] - dv0);

            double dt0 = k[iv][jdt][idt1][iidt1+i];
            dt = dt0 + rdt1*(k[iv][jdt][jdt1][jjdt1+i] - dt0);

            tmpK[i] = base + rv*(dv-base) + rdt*(dt-base);
        }
    }

    for (size i=0;i<tl;i++) {
        base = PSP0[iv][idt][fi][idt0+i];
        tmpPSP0[i] = base + rv*(PSP0[jv][idt][fi][iidt+i]-base) + rdt*(PSP0[iv][jdt][fi][jjdt+i]-base) + rf*(PSP0[iv][idt][fj][iidt+i]-base);
    }

    if (tshift!=1) {
        for (size i=1;i<tl;i++) {
            base = tmpK[i-1]*tmpPSP0[i-1]*PSP1[i-1];
            v[vs+i] = v[vs+i] + base + tshift*(tmpK[i]*tmpPSP0[i]*PSP1[i]-base);
        }
    } else {
        for (size i=0;i<tl;i++) {
            v[vs+i] = v[vs+i] + tmpK[i]*tmpPSP0[i]*PSP1[i];
        }
    }
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

unsigned int bilinear_HH(std::vector<double> &v, std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &hE, std::vector<double> &hI, std::vector<double> &m, std::vector<double> &n, std::vector<double> &h, std::vector<size> &cross, NeuroLib neuroLib, Neuron neuron, double run_t, double ignore_t, double &tsp, double pairs[], double tau_er, double tau_ed, double tau_ir, double tau_id, double vCross, double vBack, double tref, int rk, int afterSpikeBehavior, bool spikeShape, bool kVStyle){
    std::string kType1, kType2, kType;
    double f1, f2, vTarget, dtTarget, tshift, shift;
    size ts, tl, vs, ve, i, j, k, i_b = 0;
    size nE = neuroLib.nE, nI = neuroLib.nI, nv = neuroLib.nv, ndt = neuroLib.ndt;
    double tstep = neuroLib.tstep;
    size t1 = neuroLib.l0;
    size nt0 = neuroLib.nt0;
    size l0 = nt0 - 1;
    size t0 = l0 * tstep;

    size nt = neuroLib.idtRange[ndt-1]+1;
    size l1 = nt - 1;
    size t1 = l1 * tstep;

    size tb = t1 - ignore_t;
    size lb = static_cast<size>(tb/tstep);
    size nb = lb + 1;

    size run_lt = static_cast<size>(run_t/step);
    size run_nt = run_lt + 1;

    bool crossed;
    size ncross = 0;
    unsigned int spikeCount = 0;

    double *tmp = new double[nt0];
    double *tmpK = new double[nt0];
    double *tmpPSP = new double[nt0];
    double *currentPSP = new double[nt0];

    double vinit = v[0];
    v.assign(runt_nt, neuron.vReset);
    if (l0 > run_lt) {
        tl = run_nt;
    } else {
        tl = nt0;
    }
    interpVinit(v,0,neuroLib.vLeak,neuroLib.vRange,nv,vinit,1,0,tl,tmp);
    if (neuron.tin.size() > 0) {
        vs = static_cast<size>(neuron.tin[0]/tstep);
        for (i=0; i<neuron.tin.size(); i++) {
            // linear
            if (neuron.tin[i] > run_t) {
                break;
            }
            f2 = neuron.fStrength[neuron.inID[i]];
            vTarget = interpTar(v,vs,tstep,neuron.tin[i],tshift);
            shift = 1-tshift/tstep;

            if (vs+l0 > run_lt) {
                tl = run_nt - vs;
            } else {
                tl = nt0
            }

            if (f2 > 0) {
                interpPSP0(v,vs, neuroLib.EPSP0, neuroLib.vRange, neuroLib.fE, nv, nE, vTarget, f2, shift,0,tl,currentPSP);
                kType2 = "1";
            } else  {
                interpPSP0(v,vs, neuroLib.IPSP0, neuroLib.vRange, neuroLib.fI, nv, nI, vTarget,-f2, shift,0,tl,currentPSP);
                kType2 = "0";
            }

            if (i<neuron.tin.size()-1) {
                ve = static_cast<size>(neuron.tin[i+1]/tstep);
                if (ve + 1 > vs + nt) {
                    vc = vs + nt;
                } else {
                    vc = ve + 1;
                }
            } else {
                vc = run_nt;
            }

            for (k=i; k>i_b; k--) {
                j = k-1; // prevent negative for unsigned int iteration number. 
                dtTarget = neuron.tin[i]-neuron.tin[j];
                if (dtTarget > tb) break;
                dtTarget = dtTarget/tstep; // for interp along dtRange
                if (neuron.tin[j] + tb > run_t) {
                    tl = run_nt-vs;
                } else {
                    tl = static_cast<size>(ceil(nb - dtTarget));
                }
                
                if (!kVStyle) {
                    f1 = neuron.fStrength[neuron.inID[j]];
                    if (f1 > 0) kType1 = "1";
                    else kType1 = "0";
                    kType = kType1 + kType2;
                    switch (std::stoi(kType,nullptr,10)) {
                        case 11: //EE
                            interpK(v, vs, neuroLib.kVEE, neuroLib.EPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fE, nv, ndt, nE, vTarget, dtTarget, dtTarget, f1, shift, tl, tmpK, tmpPSP);
                            break;
                        case 10: //EI
                            interpK(v, vs, neuroLib.kVEI, neuroLib.EPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fE, nv, ndt, nI, vTarget, dtTarget, dtTarget, f1, shift, tl, tmpK, tmpPSP);
                            break;
                        case 1:  //IE
                            interpK(v, vs, neuroLib.kVIE, neuroLib.IPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fI, nv, ndt, nE, vTarget, dtTarget, dtTarget, -f1, shift, tl, tmpK, tmpPSP);
                            break;
                        case 0:  //II
                            interpK(v, vs, neuroLib.kVII, neuroLib.IPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fI, nv, ndt, nI, vTarget, dtTarget, dtTarget, -f1, shift, tl, tmpK, tmpPSP);
                            break;
                    }
                } else {
                    interpkV(v, vs, neuroLib.kV, neuroLib.vRange, neuroLib.idtRange, nv, ndt, vTarget, dtTarget, dtTarget, shift, tl, tmpK);
                    
                }
            } 
            //check spike 
            crossed = false;
            ith = i;
            for (k=vs;k<vc;k++) {
                while (v[k] > vCross) {
                    crossed = true;
                    ncross = ncross +1;
                    // come in and out in multiples of tstep 
                    std::cout << "crossed at " << k*tstep << ", start ith " << ith << std::endl;
                    if (spikeShape) {
                        spiked = nearThretyleshold(neuron, neuroLib, v, gE, gI, hE, hI, m, n, h, k, nt+1, tstep, pairs, tau_er, tau_ed, tau_ir, tau_id, tsp, v[k], vs, ith, rk, vBack);
                        cross.push_back(gE.size());
                    } else {
                        size ith_old = ith; 
                        size itref = static_cast<size>(round(neuron.tref/tstep));
                        spiked = 1;
                        double tmpTsp = k*neuroLib.tstep+neuron.tref/2;
                        if (tmpTsp <= run_t) {
                            tsp.push_back(tmpTsp);
                        }
                        vs = k + itref;
                        for (j=0;j<itref;j++) {
                            if (k+j < run_nt) {
                                v[k+j] = vCross;
                            } else {
                                vs = k+j-1;
                                break;
                            }
                        }
                        v[vs] = neuron.vRest;
                        for (ith = ith_old; ith < nin; ith++) {
                            if (neuron.tin[ith] > vs*tstep) {
                                break;                
                            }
                        }
                        ith--;
                        cross.push_back(itref);
                    }
                    std::cout << "backed at " << vs*tstep << ", end ith " << ith << std::endl;
                    std::cout << "v[k] " << v[vs]  << " vCross " << vCross << " vBack " << vBack << std::endl;
                    if (spiked){
                        spikeCount = spikeCount + 1;
                        neuron.tsp.push_back(tsp);
                    } 
                    if (vs + nt0 <= nt) {
                        tl = nt0;
                    } else {
                        tl = nt-vs+1;
                    }

                    for (j=1; j<tl;j++) {
                        v[vs+j] = neuron.vRest;
                    }

                    interpVinit(v,vs,neuroLib.vLeak,neuroLib.vRange,nv,v[vs],1,0,tl,tmp);
                    if (afterSpikeBehavior) {
                        vTarget = interpTar(v,vs,tstep,neuron.tin[j],tshift);
                        shift = 1-tshift/tstep;

                        for (i=ith+1;i>0; i--) {
                            j = i-1;
                            it = round(neuron.tin[j]/tstep);
                            dtTarget = vs - it;
                            if ( l1<= dtTarget) {
                                break;
                            }
                            if (run_t < neuron.tin[j] + t1)
                                tl = run_nt-vs;
                            else tl = static_cast<size>(ceil(nt - dtTarget))+1;
                            f2 = neuron.fStrength[neuron.inID[j]];
                            if (f2 > 0) {
                                interpPSP(v,vs, neuroLib.EPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fE, nv, ndt, nE, vTarget, dtTarget, f, 1,0,tl,currentPSP);
                                kType2 = "1";
                            } else {
                                interpPSP(v,vs, neuroLib.IPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fI, nv, ndt, nI, vTarget, dtTarget,-f, 1,0,tl,currentPSP);
                                kType2 = "0";
                            }
                            if (afterSpikeBehavior == 2) {
                                i_b = 0;
                                cout << " adding afterspike bPSP " << endl;
                                for (k=i; k>0; k--) {
                                    j = k-1; // prevent negative for unsigned int iteration number. 
                                    jt = neuron.tin[j]/tstep;
                                    dtTarget = it - jt;
                                    dtTarget1 = vs - jt;
                                    if (dtTarget > lb) {
                                        break;
                                    } 
                                    if (neuron.tin[j] + tb > run_t) {
                                        tl = run_nt - vs;
                                    } else { 
                                        tl = static_cast<size>(round(nb - dtTarget1));
                                    }
                                    //cout << "vs + tl " << vs << " + " << tl << " < " << run_nt << endl;
                                    if (!kVStyle) {
                                        f1 = neuron.fStrength[neuron.inID[j]];
                                        if (f1 > 0) kType1 = "1";
                                        else kType1 = "0";
                                        kType = kType1 + kType2;
                                        switch (std::stoi(kType,nullptr,10)) {
                                            case 11: //EE
                                                interpK(v, vs, neuroLib.kVEE, neuroLib.EPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fE, nv, ndt, nE, vTarget, dtTarget, dtTarget1, f1, shift, tl, tmpK, tmpPSP);
                                                break;
                                            case 10: //EI
                                                interpK(v, vs, neuroLib.kVEI, neuroLib.EPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fE, nv, ndt, nI, vTarget, dtTarget, dtTarget1, f1, shift, tl, tmpK, tmpPSP);
                                                break;
                                            case 1:  //IE
                                                interpK(v, vs, neuroLib.kVIE, neuroLib.IPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fI, nv, ndt, nE, vTarget, dtTarget, dtTarget1, -f1, shift, tl, tmpK, tmpPSP);
                                                break;
                                            case 0:  //II
                                                interpK(v, vs, neuroLib.kVII, neuroLib.IPSP, currentPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fI, nv, ndt, nI, vTarget, dtTarget, dtTarget1, -f1, shift, tl, tmpK, tmpPSP);
                                                break;
                                        }

                                    } else {
                                        interpkV(v, vs, neuroLib.kV, neuroLib.vRange, neuroLib.idtRange, nv, ndt, vTarget, dtTarget, dtTarget1, shift, tl, tmpK);
                                    }
                                } 
                            }
                        }
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            if (ve > vs + l1){
                                vc = vs + nt; 
                            } else {
                                vc = ve + 1;
                            }
                        } else {
                            vc = run_nt;
                        }
                        k = vc-1;
                        for (k=vs; k<vc; k++) {
                            if (v[k] > vCross) {
                                k = j;
                                std::cout << " need to cross at " << k << std::endl;
                                break;
                            }
                        }
                    } else {
                        cout << " no readjust input, only leakage" << endl;
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                        }
                        break;
                    }
                }
                if (crossed) {
                    i = ith;
                    cout << " last input during cross " << i << endl;
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

unsigned int linear_HH(std::vector<double> &v, std::vector<double> &gE, std::vector<double> &gI, std::vector<double> &hE, std::vector<double> &hI, std::vector<double> &m, std::vector<double> &n, std::vector<double> &h, std::vector<size> &cross, NeuroLib neuroLib, Neuron neuron, double run_t, double ignore_t, double &tsp, double pairs[], double tau_er, double tau_ed, double tau_ir, double tau_id, double vCross, double vBack, double tref, int rk, int afterSpikeBehavior, bool spikeShape){
    double f, vTarget, dtTarget, tshift, shift;
    size ts, tl, vs, ve, i, j, k, i_b = 0;
    size nE = neuroLib.nE, nI = neuroLib.nI, nv = neuroLib.nv, ndt = neuroLib.ndt;
    double tstep = neuroLib.tstep;
    size t1 = neuroLib.l0;
    size nt0 = neuroLib.nt0;
    size l0 = nt0 - 1;
    size t0 = l0 * tstep;

    size nt = neuroLib.idtRange[ndt-1]+1;
    size l1 = nt - 1;
    size t1 = l1 * tstep;

    size tb = t1 - ignore_t;
    size lb = static_cast<size>(tb/tstep);
    size nb = lb + 1;

    size run_lt = static_cast<size>(run_t/step);
    size run_nt = run_lt + 1;

    bool crossed;
    size ncross = 0;
    unsigned int spikeCount = 0;

    double *tmp = new double[nt0];
    double *tmpK = new double[nt0];
    double *tmpPSP = new double[nt0];
    double *currentPSP = new double[nt0];

    double vinit = v[0];
    v.assign(runt_nt, neuron.vReset);
    if (l0 > run_lt) {
        tl = run_nt;
    } else {
        tl = nt0;
    }
    interpVinit(v,0,neuroLib.vLeak,neuroLib.vRange,nv,vinit,1,0,tl,tmp);
    if (neuron.tin.size() > 0) {
        vs = static_cast<size>(neuron.tin[0]/tstep);
        for (i=0; i<neuron.tin.size(); i++) {
            // linear
            if (neuron.tin[i] > run_t) {
                break;
            }
            f = neuron.fStrength[neuron.inID[i]];
            vTarget = interpTar(v,vs,tstep,neuron.tin[i],tshift);
            shift = 1-tshift/tstep;

            if (vs+l0 > run_lt) {
                tl = run_nt - vs;
            } else {
                tl = nt0
            }

            if (f > 0) {
                interpPSP0(v,vs, neuroLib.EPSP0, neuroLib.vRange, neuroLib.fE, nv, nE, vTarget, f, shift,0,tl,currentPSP);
            } else  {
                interpPSP0(v,vs, neuroLib.IPSP0, neuroLib.vRange, neuroLib.fI, nv, nI, vTarget,-f, shift,0,tl,currentPSP);
            }

            if (i<neuron.tin.size()-1) {
                ve = static_cast<size>(neuron.tin[i+1]/tstep);
                if (ve + 1 > vs + nt) {
                    vc = vs + nt;
                } else {
                    vc = ve + 1;
                }
            } else {
                vc = run_nt;
            }
            //check spike 
            crossed = false;
            ith = i;
            for (k=vs;k<vc;k++) {
                while (v[k] > vCross) {
                    crossed = true;
                    ncross = ncross +1;
                    // come in and out in multiples of tstep 
                    std::cout << "crossed at " << k*tstep << ", start ith " << ith << std::endl;
                    if (spikeShape) {
                        spiked = nearThretyleshold(neuron, neuroLib, v, gE, gI, hE, hI, m, n, h, k, nt+1, tstep, pairs, tau_er, tau_ed, tau_ir, tau_id, tsp, v[k], vs, ith, rk, vBack);
                        cross.push_back(gE.size());
                    } else {
                        size ith_old = ith; 
                        size itref = static_cast<size>(round(neuron.tref/tstep));
                        spiked = 1;
                        double tmpTsp = k*neuroLib.tstep+neuron.tref/2;
                        if (tmpTsp <= run_t) {
                            tsp.push_back(tmpTsp);
                        }
                        vs = k + itref;
                        for (j=0;j<itref;j++) {
                            if (k+j < run_nt) {
                                v[k+j] = vCross;
                            } else {
                                vs = k+j-1;
                                break;
                            }
                        }
                        v[vs] = neuron.vRest;
                        for (ith = ith_old; ith < nin; ith++) {
                            if (neuron.tin[ith] > vs*tstep) {
                                break;                
                            }
                        }
                        ith--;
                        cross.push_back(itref);
                    }
                    std::cout << "backed at " << vs*tstep << ", end ith " << ith << std::endl;
                    std::cout << "v[k] " << v[vs]  << " vCross " << vCross << " vBack " << vBack << std::endl;
                    if (spiked){
                        spikeCount = spikeCount + 1;
                        neuron.tsp.push_back(tsp);
                    } 
                    if (vs + nt0 <= nt) {
                        tl = nt0;
                    } else {
                        tl = nt-vs+1;
                    }

                    for (j=1; j<tl;j++) {
                        v[vs+j] = neuron.vRest;
                    }

                    interpVinit(v,vs,neuroLib.vLeak,neuroLib.vRange,nv,v[vs],1,0,tl,tmp);
                    if (afterSpikeBehavior) {
                        vTarget = interpTar(v,vs,tstep,neuron.tin[j],tshift);
                        shift = 1-tshift/tstep;

                        for (i=ith+1;i>0; i--) {
                            j = i-1;
                            it = round(neuron.tin[j]/tstep);
                            dtTarget = vs - it;
                            if ( l1<= dtTarget) {
                                break;
                            }
                            if (run_t < neuron.tin[j] + t1)
                                tl = run_nt-vs;
                            else tl = static_cast<size>(ceil(nt - dtTarget))+1;
                            f = neuron.fStrength[neuron.inID[j]];
                            if (f > 0) {
                                interpPSP(v,vs, neuroLib.EPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fE, nv, ndt, nE, vTarget, dtTarget, f, 1,0,tl,currentPSP);
                            } else {
                                interpPSP(v,vs, neuroLib.IPSP, neuroLib.vRange, neuroLib.idtRange, neuroLib.fI, nv, ndt, nI, vTarget, dtTarget,-f, 1,0,tl,currentPSP);
                            }
                        }
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                            if (ve > vs + l1){
                                vc = vs + nt; 
                            } else {
                                vc = ve + 1;
                            }
                        } else {
                            vc = run_nt;
                        }
                        k = vc-1;
                        for (k=vs; k<vc; k++) {
                            if (v[k] > vCross) {
                                k = j;
                                std::cout << " need to cross at " << k << std::endl;
                                break;
                            }
                        }
                    } else {
                        cout << " no readjust input, only leakage" << endl;
                        if (ith<neuron.tin.size()-1) {
                            ve = static_cast<size>(neuron.tin[ith+1]/tstep);
                        }
                        break;
                    }
                }
                if (crossed) {
                    i = ith;
                    cout << " last input during cross " << i << endl;
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
