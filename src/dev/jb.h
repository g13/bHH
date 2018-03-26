#ifndef JB_H 
#define JB_H
#include <cmath>
#include <vector>
#include "typedefs.h"
#include "nStructs.h"
#include "NeuronStruct.h"

using std::vector;
using std::cout;
using std::endl;
namespace jb{
    const bool debug = false;
    const bool debug2 = false;
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
}

inline void add_relation_between_new_and_ith_input(Input &input, size inew, size i, double *idtRange) {
    double t = input.t[inew] - input.t[i];
    if (t < 0) {
        std::cout << i << ", " << inew << std::endl;
        std::cout << input.t[i] << ", " << input.t[inew] << std::endl;
    }
    double r_;
    size i_,j_;
    //size idt = static_cast<size>(t);
    jb::getNear(idtRange, neuroLib.ndt, t, r_,i_,j_);
    input.bir[i].dTijr.push_back(IJR(i_,j_,r_));
    //input.bir[i].Tijr.push_back(IJR(idt,idt+1,t-idt));
}

inline double linear_interp_tMax(double ***tMax, IJR &v, IJR &dT, size i) {
    double base = tMax[v.i][dT.i][i];
    return round(base + v.r *  (tMax[v.j][dT.i][i] -base)
                      + dT.r * (tMax[v.i][dT.j][i] -base));
}

inline void add_new_input_info(double *vRange, size nv, Input &input, Cross &cross, double t0, double ***tMax, double v, size corrSize, size ID) {
    double t = input.t.back();
    double r_;
    size i_, j_;
    jb::getNear(vRange, nv, v, r_,i_,j_);
    input.Vijr.push_back(IJR(i_,j_,r_));
    input.dTijr.push_back(IJR(0,1,0));
    input.dt.push_back(t);
    input.Tmpijr.push_back(IJR(0,1,0));
    input.cCross.push_back(0);
    input.tMax.push_back(linear_interp_tMax(tMax, input.Vijr.back(), input.dTijr.back(),ID));
    input.bir.push_back(BilinearRelationships(corrSize));
}

inline double add_vinit_contribution(double **vLeak, IJR &v, double dt) { 
    size i = static_cast<size>(round(dt));
    double base = vLeak[v.i][i];
    return base + v.r * (vLeak[v.j][i] -base);
}

inline void move_corr_window(vector<double> &input, size &tail, double t_head, double tol_tb, double tstep) {
    double dt;
    dt = t_head - round(input[tail]/tstep);
    while (dt > tol_tb) {
        tail++;
        dt = t_head - round(input[tail]/tstep);
    }
}

inline void reverse_corr_window(vector<double> &input, size &tail, double t_head, double tol_tb, double tstep) {
    double dt;
    do {
        tail--;
        if (tail<0) {
            break;
        }
        dt = t_head - round(input[tail]/tstep);
    } while (dt <= tol_tb);
    tail++;
}

inline double linear_interp_kV(double ******kV, IJR &v, IJR &ijdT, size dT, size it, size *idtRange, size ndt, size i, size j, bool debug) {
    if (dT == 0) {
        size jt = idtRange[ijdT.j] + it;
        it = idtRange[ijdT.i] + it;
        if (debug) {
            cout << "     dT " << dT << ": " << ijdT.i << ", " << ijdT.j << ", r" << ijdT.r << endl;
            cout << "     it " << it <<", jt" << jt << endl;
        }
        double base = kV[v.i][ijdT.i][i][j][ijdT.i][it];
        double dv = kV[v.j][ijdT.i][i][j][ijdT.i][it] - base;
        double dt = kV[v.i][ijdT.j][i][j][ijdT.j][jt] - base;
        if (debug) {
            cout << "     base " << base << ", dv " << dv << ", dt " << dt << endl;
        }
        return base + v.r * dv + ijdT.r * dt;
    } else {
        size idt0, jdt0, idt1, jdt1;
        double rdt0, rdt1;
        jb::getNear(idtRange,ndt,idtRange[ijdT.i] + dT, rdt0, idt0, jdt0);
        jb::getNear(idtRange,ndt,idtRange[ijdT.j] + dT, rdt1, idt1, jdt1);
        size iidt0 = idtRange[idt0]+it;
        size jjdt0 = idtRange[jdt0]+it;
        size iidt1 = idtRange[idt1]+it;
        size jjdt1 = idtRange[jdt1]+it;
        if (debug) {
            cout << "     dT " << dT << ": " << iidt0 <<", " << jjdt0 << ", r" << rdt0 << "; " << iidt1 << ", " << jjdt1 << ", r" << rdt1 << " it " << it << endl;
        }

        double base0 = kV[v.i][ijdT.i][i][j][idt0][iidt0];
        double base = base0 + rdt0*(kV[v.i][ijdT.i][i][j][jdt0][jjdt0]-base0);

        double dv0 = kV[v.j][ijdT.i][i][j][idt0][iidt0];
        double dv = dv0 + rdt0*(kV[v.j][ijdT.i][i][j][jdt0][jjdt0]-dv0);

        double dt0 = kV[v.i][ijdT.j][i][j][idt1][iidt1];
        double dt = dt0 + rdt1*(kV[v.i][ijdT.j][i][j][jdt1][jjdt1]-dt0);

        if (debug) {
            cout << base0 << ", " << dv0 << ", " << dt0 << endl;
            cout << "inds: "<< v.i << ", "<< ijdT.j<< ", " << i << ", " << j << ", " << idt1 << ", " << iidt1 << endl;
        }
        return base + v.r * (dv-base) + ijdT.r * (dt-base);
    }
}

inline double linear_interp_PSP(double ****PSP, IJR &v, IJR &dT, size ID, size idt, size *idtRange) {
    size iidt = idt + idtRange[dT.i];
    size jjdt = idt + idtRange[dT.j];
    //cout << " dT.i " << dT.i << endl;
    //cout << " dT.r " << dT.r << endl;
    double base = PSP[v.i][dT.i][ID][iidt];
    return base + v.r *  (PSP[v.j][dT.i][ID][iidt] -base)
                + dT.r * (PSP[v.i][dT.j][ID][jjdt] -base);
}

inline void add_input_i_contribution(size i, size idt, NeuroLib &neuroLib, Input &input, double &v) {
    //cout << "synapse " << input.ID[i] << " contributing " << input.t[i] << endl;
    if (input.ID[i] < neuroLib.nE) {
        v = v + linear_interp_PSP(neuroLib.sEPSP, 
                                  input.Vijr[i], input.dTijr[i], 
                                  input.ID[i], idt, neuroLib.idtRange);
    else {
        v = v + linear_interp_PSP(neuroLib.sIPSP, 
                                  input.Vijr[i], input.dTijr[i], 
                                  input.ID[i]-nE, idt, neuroLib.idtRange);
    }
}
inline void add_input_i_j_bilinear_contribution(Input &input, NeuroLib &neuroLib, size i, size j, size it, double &v, bool debug) {
    size k = j-i-1;
    if (input.cCross[j] > input.cCross[i]) {
        cout << "i.cCross " << input.cCross[i] << " < j.cCross " << input.cCross[j] << endl;
        cout << i << " t " << input.t[i] << j << " t " << input.t[j] <<  "-> " << input.t[j] - input.t[i] << endl;
        assert(input.cCross[j] <= input.cCross[i]);
    }
    v += linear_interp_kV(neuroLib.kV, input.Vijr[j], input.bir[i].dTijr[k], 
                          input.dt[j]-input.t[j], it, neuroLib.idtRange, neuroLib.ndt, input.ID[i], input.ID[j], debug);
}

inline double find_v_at_t(Input &input, NeuroLib &neuroLib, Cross &cross, size head, size tail_l, size tail_b, double t, double tCross, double tol_tl, double v, bool debug) {
    long i;
    size j, it;
    double dt, dt0;
    dt = t - tCross;
    // v = neuron.vReset unless
    if (dt < tol_tl) {
        v = add_vinit_contribution(neuroLib.vLeak, cross.vCross.back(), dt);
    } 
    for (i=head; i>=tail_l; i--) {
        dt0 = t - input.t[i];
        if (dt0 > tol_tl) break;
        dt = t - input.dt[i];
        it = static_cast<size>(round(dt));
        add_input_i_contribution(i,it,neuroLib,input,v);
        //cout << " i " << i << " cCross" << input.cCross[i] << endl;
        if (i < tail_b) continue;
        for (j=head; j>i; j--) {
            //cout << "   j " << j << " cCross" << input.cCross[j] << endl;
            dt = t - input.dt[j];
            it = static_cast<size>(round(dt));
            add_input_i_j_bilinear_contribution(input, neuroLib, i, j, it, v, debug);
        }
    }
    return v;
}
inline double parabola(double t_left, double v_left, double t_right, double v_right, double t_mid, double v_mid, double vThres) {
    double t1 = t_mid-t_left;
    double t1_2 = t1*t1;
    double t2 = t_right-t_left;
    double t2_2 = t2*t2; 
    double denorm = t1_2*t2-t2_2*t1;
    assert(denorm != 0.0);
    v_right = v_right-v_left;
    v_mid = v_right-v_left;
    vThres = vThres-v_left;
    double a = (v_mid*t2-v_right*t1)/denorm;
    double b = (v_right*t1_2-v_mid*t2_2)/denorm;
    return (-b+sqrt(b*b+4*a*vThres))/(2*a) + t_left;
}
inline void getLR(double &v_left, double v1, double v2, double &v_right, double &t_left, double t_cross1, double t_cross2, double &t_right, double vTar) {
    if (vTar < v1) {
        v_right = v1;
        t_right = t_cross1;
    } else {
        if (vTar < v2) {
            v_left = v1;
            t_left = t_cross1;
            v_right = v2;
            t_right = t_cross2;
        } else {
            v_left = v2;
            t_left = t_cross2;
        }
    }
}

double interp_for_t_cross(double v_right, double v_left, double t_right, double t_left, size head, size tail_l, size tail_b, double tCross, double tol_tl, NeuroLib &neuroLib, Cross &cross, Input &input, double v0, double vtol, double vC, double &v, bool debug) {
    size ival = 0;
    double v1, v2
    double t_cross1, t_cross2, t_cross;
    do {
        assert(vC >= v_left);
        if (abs(t_right - t_left)<1e-14) {
            v = v_right;
            return t_right;
        }
        assert(t_right > t_left);
        t_cross1 = t_left + (vC - v_left)/(v_right-v_left)*(t_right-t_left);
        //t_cross1 = ceil(t_cross1);
        assert(t_cross1 >= t_left);
        if (jb::debug2) {
            std::cout << " find v at " << t_cross1 << std::endl;
            std::cout <<"t: " << t_left << ", " << t_right << std::endl;
        }
        v1 = find_v_at_t(input, neuroLib, cross, head, tail_l, tail_b, t_cross1, tCross, tol_tl, v0, debug);
        if (jb::debug2) {
            std::cout << v_left <<  ", " << v1 << ", " << v_right << std::endl;
        }
        ival++;
        if (fabs(v1-vC)<vtol) {
            v = v1;
            t_cross = t_cross1;
            break;
        }
        if (t_cross1-t_left <= 1 || t_right-t_cross1 <=1) {
            if (t_cross1-t_left <= 1) {
                v = v1; 
            } else {
                v = v_right;
            }
            t_cross = t_cross1;
            if (jb::debug2) {
                std::cout << "temp resol reached " << std::endl;
            }
            break;
        }
        // solve for a, b of f(t)-v_left = a(t-t_left)^2 + b(t-t_left);
        // solve for t when a(t-t_left)^2 + b(t-t_left) = (vThres - v_left)
        t_cross2 = parabola(t_left,v_left,t_right,v_right,t_cross1, v1, vC);
        //t_cross2 = ceil(t_cross2);
        if (jb::debug2) {
            std::cout << " v: " << v_left << ", " << v1 << ", " << v_right << std::endl;
            std::cout << " t: " << t_left << ", " << t_cross1 << ", " << t_right << std::endl;
            std::cout << "t_cross2 " << t_cross2 << std::endl;
            assert(t_cross2 >= t_left);
            assert(t_cross2 <= t_right);
        }
        // if t interval smaller than sample temp resolution, apply the quadratic interpolation as solution
        if (jb::debug2) {
            std::cout << " find v at " << t_cross2 << std::endl;
        }
        v2 = find_v_at_t(input, neuroLib, cross, head, tail_l, tail_b, t_cross2, tCross, tol_tl, v0, debug);
        //if (jb::debug2) {
            ival++;
        //}
        if (fabs(v2-vC)<vtol) {
            v = v2;
            t_cross = t_cross2;
            break;
        }
        if (jb::debug2) {
            std::cout << " v: " << v_left << ", " << v1 << ", " << v2 << ", " << v_right << std::endl;
            std::cout << " t: " << t_left << ", " << t_cross1 << ", " << t_cross2 << ", " << t_right << std::endl;
        }
        if (v2>v1) {
            getLR(v_left,v1,v2,v_right,t_left,t_cross1,t_cross2,t_right, vC);
        } else {
            getLR(v_left,v2,v1,v_right,t_left,t_cross2,t_cross1,t_right, vC);
        } 
        if (jb::debug2) {
            std::cout << " v: " << v_left << ", " << v_right << std::endl;
            std::cout << " t: " << t_left << ", " << t_right << std::endl;
        }
    } while (true);
    //if (jb::debug2) {
        std::cout << ival << " evaluations" << std::endl;
    //}
    return t_cross;
}

bool check_crossing(Input &input, NeuroLib &neuroLib, Cross &cross, Neuron &neuron, double tol_tl, double tol_tb, double end_t, size tail_l, size tail_b, size head, jND &jnd, double &t_cross, double vC, bool debug) {
    // check all tmax after head and find the upper limit of vmax
    size i, j, it;
    double t, dt;
    size nt = end_t/neuroLib.tstep;
    double v, v_pass; 
    double vmax = jnd.v.back();
    double tmax = jnd.t.back();
    size vi_tail;
    size vi_tail_max;
    double tCross = cross.tCross.back();
    for (i=tail_l; i<=head; i++) {
        t = input.tMax[i]+input.t[i];
        vi_tail = 0;
        if (t <= jnd.t.back() || input.ID[i] >= neuroLib.nE) {
            // ignore inh input and input that uncorred after cross
            continue;
        }
        // ignore tmax that go beyond next input
        if (head < neuron.tin.size()-1) {
            if (t >= neuron.tin[head+1]/neuroLib.tstep) {
                continue;
            }
        } else {
            if ( t > nt ){
                t = nt;
            }
        }
        // initialize with leak from last cross
        dt = t-tCross;
        if (dt < tol_tl) {
            v = add_vinit_contribution(neuroLib.vLeak, cross.vCross.back(), dt);
        } else {
            v = neuron.vReset;
        }
        for (j=tail_l; j<=head; j++) {
            dt = t - input.t[j];
            if (dt > tol_tl) {
                continue;
            } else { 
                // get tail for tmax
                if ( dt < tol_tb ) {
                    if (!vi_tail) {
                        vi_tail = j;
                    }
                }
            }
            dt = t - input.dt[j];
            it = static_cast<size>(round(dt));
            add_input_i_contribution(j, it, neuroLib, input, v);
        }
        if (v > vmax) {
            vmax = v;
            tmax = t; 
            vi_tail_max = vi_tail;
        }
    }
    if (vmax > vCross) {
        if (jb::debug2) {
            std::cout << "linear vmax " << vmax <<  " > vThres" << std::endl;
        }
        // perform linear interp iteration until tolerance reaches
        if (jb::debug) {
            std::cout << "l+b vmax > vThres, find t and v for cross" << std::endl;
        }
        if (fabs(vmax-vC)>vtol && tmax-jnd.t.back()>1) {
            t_cross = interp_for_t_cross(vmax, jnd.v.back(), tmax, jnd.t.back(), head, tail_l, tail_b, tCross, tol_tl, neuroLib, cross, input, neuron.vReset, neuron.vTol, vC, v_pass, debug);
            assert(t_cross > jnd.t.back());
            jnd.t.push_back(t_cross);
            jnd.v.push_back(v_pass);
        } else {
            t_cross = tmax;
            jnd.t.push_back(tmax);    
            jnd.v.push_back(vmax);    
        }
        return true;
    }
    return false;
}

bool update_vinit_of_new_input_check_crossing(Input &input, Cross &cross, NeuroLib &neuroLib, Neuron &neuron, size head, size tail_l, size tail_b, double tol_tl, double tol_tb, double end_t, jND &jnd, double &t_cross, double vC, size corrSize, bool debug) {
    long i;
    size j;
    double dt;
    double tstep = neuroLib.tstep;
    double v_pass;
    double ir;
    size idt;
    double tCross = cross.tCross.back();
    double v;
    //if (jb::debug) {
    //    cout << " adding " << head << " input t " << input.t[head] << " synapse ID " << neuron.inID[head] << endl;
    //}
    input.ID.push_back(neuron.inID[head]);
    dt = input.t[head] - tCross;
    if (jb::debug) {
        if (dt < -1e-14 ) {
            std::cout << dt << " < 0 " << std::endl;
            std::cout << head << std::endl;
            std::cout << neuron.tin[head] << std::endl;
            std::cout << input.t[head] << std::endl;
            std::cout << tCross*tstep << std::endl;
            assert(dt>=0);
        }
    }
    // initialize to 
    if (dt < tol_tl) {
        v = add_vinit_contribution(neuroLib.vLeak, cross.vCross.back(), dt);
    } else {
        v = neuron.vReset;
    }
    if (debug) {
        cout << " last v " << jnd.v.back() << endl;
        if (head-1 < tail_l && tail_l != 0) {
            cout << "first v after cross ";
        } else {
            cout << " start with leak";
        }
        cout << v << endl;
    }
    //cout << " tailing input " << tail_l << endl;
    double vtmp = v;
    double dvtmp;
    for (i=head-1; i>=tail_l; i--) {
        if (debug) {
            dvtmp = 0;
        }
        dt = input.t[head] - input.dt[i];
        input.bir[i].idt.push_back(static_cast<size>(round(dt)));
        input.bir[i].ID.push_back(neuron.inID[head]);
        add_input_i_contribution(i,input.bir[i].idt.back(),neuroLib,input,v);
        if (debug) {
            cout << " " << i << " input contributing " << v-vtmp << ", synapse ID " << input.ID[i] << endl;
            vtmp = v;
        }
        if (i<tail_b) continue;
        for (j=head-1; j>i; j--) {
            ir = head-j-1;
            add_input_i_j_bilinear_contribution(input, neuroLib, i, j, input.bir[j].idt[ir], v, debug);
            if (debug) {
                dvtmp += v-vtmp;
                cout << "    +" << j << " bilinear contributing " << v-vtmp << ", synapse ID " << input.ID[j] << " with dti " << input.bir[i].dTijr[j-i-1].i << ", r" << input.bir[i].dTijr[j-i-1].r << " at vi " << input.Vijr[j].i << ", r" << input.Vijr[j].r << endl;
                vtmp = v;
            }
        }
        if (debug) {
            cout << " bilinear total: " << dvtmp << endl;
            cout << " v: " << v << endl;
        }
        
    }
    // update head's v and tmax
    jnd.t.push_back(input.t[head]);
    jnd.v.push_back(v);
    //if (jb::debug) {
    //    cout << "v_" << head << " = " << jnd.v.back() << endl;
    //}
    // check if vmax at left bound very unlikely
    if (v > vC) {
        if(jb::debug2) {
            std::cout << " crossing upon input not very probable at low input rate " << std::endl;
        }
        jnd.t.pop_back();
        jnd.v.pop_back();
        // perform linear interp iteration until tolerance reaches
        cout << "t: " << jnd.t.back()  << " < " << input.t[head] << endl;
        t_cross = interp_for_t_cross(v, jnd.v.back(), input.t[head], jnd.t.back(), head-1, tail_l, tail_b, tCross, tol_tl, neuroLib, cross, input, neuron.vReset, neuron.vTol, vC, v_pass, debug);
        if (t_cross-input.t[head] > 1e-10) {
        
            assert(t_cross<=input.t[head]);
        }
        jnd.t.push_back(t_cross);
        jnd.v.push_back(v_pass);
        input.t.pop_back();
        input.ID.pop_back();
        for (i=head-1; i>=tail_l; i--) {
            input.bir[i].idt.pop_back();
            input.bir[i].ID.pop_back();
        }
        if (input.t.size() != head) {
            cout << "before check " << input.t.size() << " != " << head << endl;
            assert(input.t.size() == head);
        }
        return true;
    } else {
        add_new_input_info(neuroLib.vRange, neuroLib.nv, input, cross, neuron.tin[head], neuroLib.tMax, v, corrSize, neuron.inID[head]);
        for (i=head-1; i>=tail_l; i--) {
            add_relation_between_new_and_ith_input(input, head, i, neuroLib.idtRange);
        }
        if (input.t.size() != head+1) {
            cout << "before check " << input.t.size() << " != " << head+1 << endl;
            assert(input.t.size() == head+1);
        }
    }
    // check all tmax after head and find the upper limit of vmax
    return check_crossing(input, neuroLib, cross, neuron, tol_tl, tol_tb, end_t, tail_l, tail_b, head, jnd, t_cross, vC, debug);
}

void update_info_after_cross(Input &input, NeuroLib &neuroLib, Cross &cross, Neuron &neuron, double tCross, double vCross, size i_prior, size tail_l, size tail_b, size head, size corrSize, int afterCrossBehavior) {
    size i_, j_;
    double r_;
    double dt, v;
    size i, j;
    size idt;
    // inputs that only matters during cross, can fill junk with those
    size i_cross = i_prior + 1;
    size i_start = i_cross;
    if (input.t.size() != i_start) {
        cout << "before update " << input.t.size() << " != " << i_start << endl;
        assert(input.t.size() == i_start);
    }
    if (!afterCrossBehavior) {
        for (i=i_cross; i<=head; i++) {
            input.junk_this();
        }
    } else {
        if (i_cross < tail_l) {
            std::cout << " crossing part lasting longer than PSP sample length, very unlikely" << std::endl;
            cout << " dead old input " << tail_l-i_cross << endl;
            for (i=i_cross; i<tail_l; i++) {
                input.junk_this();
            }
            i_start = tail_l;
        } else {
        // update for input that come before cross and lingers after cross;
            if (jb::debug) {
                cout << "cross update starting tail " << tail_l << endl;
                cout << " # lingering old inputs: " << i_cross - tail_l << endl;
            }
            for (i=tail_l; i<i_cross; i++) {
                if (jb::debug) {
                    cout << "lingering " << i << " < " << neuron.tin.size() << endl;
                }
                dt = tCross - input.t[i];
                assert(dt>0);
                jb::getNear(neuroLib.idtRange,neuroLib.ndt,
                              dt, input.dTijr[i].r, input.dTijr[i].i, input.dTijr[i].j);
                input.dt[i] = tCross;
                input.cCross[i] = cross.nCross;
                input.Vijr[i] = cross.vCross.back();
                input.tMax[i] = linear_interp_tMax(neuroLib.tMax, cross.vCross.back(), input.dTijr[i],input.ID[i]);
                //cout << "   loop ended" << endl;
            }
        }
        // update for new input during the cross that lingers after cross
        if (jb::debug) {
            cout << " # lingering new inputs: " << head-i_start+1 << endl;
            cout << input.t.size() << " == " << input.Vijr.size();
        }
        if (input.t.size() != i_start) {
            cout << input.t.size() << " != " << i_start << endl;
            assert(input.t.size() == i_start);
        }
        for (i=i_start; i<=head; i++) {
            //std::cout << "i " << i << " < " << neuron.tin.size() << " == " << neuron.inID.size() << std::endl;
            input.t.push_back(neuron.tin[i]/neuroLib.tstep);
            assert(input.t.size() == i+1);
            input.ID.push_back(neuron.inID[i]);
            if (jb::debug) {
                cout << head << " >= " << i << " > " << i_start << endl;
                cout << i << " == " << input.t.size()-1 << endl;
                cout << input.t[i] << " == " << input.t.back() << endl;
            }
            if (input.t[i] != input.t.back()) {
                std::cout << i << " == " << input.t.size() -1 << std::endl;
                std::cout << input.t[i]*neuroLib.tstep << " == " << neuron.tin[i]/neuroLib.tstep << std::endl;
                std::cout << input.t[i]*neuroLib.tstep << " == " << input.t[i-1]*neuroLib.tstep << std::endl;
                assert(input.t[i] == input.t.back());
            }
            if ( input.t[i] - tCross > 1e-14) {
                cout << "input.t " << i << " = " << input.t[i] << " < " << tCross << endl;
                assert(tCross>=input.t[i]);
            }
            input.dt.push_back(tCross);
            dt = tCross - input.t[i];
            if (dt < -1e-14) {
                cout << tCross << " - " << input.t[i] << " = " << dt << endl;
                assert(dt >= 0);
            }
            jb::getNear(neuroLib.idtRange, neuroLib.ndt,
                           dt, r_, i_, j_);
            if (r_-1>1e-14) {
                cout << tCross  << "-" << input.t[i] << endl;
                cout << "tail " << tail_l << ", " << input.t[tail_l] << endl;
                cout << i << " head " << head << "total " << neuron.tin.size() << endl;
                cout << "dt " << dt << " i_ " << i_ << " j_ " << j_ << endl;
                cout << "r_ " << r_ << endl;
                assert(r_ <= 1);
            }
            input.dTijr.push_back(IJR(i_,j_,r_));
            input.cCross.push_back(cross.nCross);
            input.Vijr.push_back(cross.vCross.back());
            input.tMax.push_back(linear_interp_tMax(neuroLib.tMax, cross.vCross.back(), input.dTijr[i],input.ID[i]));
            input.Tmpijr.push_back(IJR(0,1,0)); 
            input.bir.push_back(BilinearRelationships(corrSize));
            if (afterCrossBehavior == 2) {
                for (j=tail_b; j<i; j++) {
                    dt = input.t[i] - input.dt[j];
                    input.bir[j].idt.push_back(static_cast<size>(round(dt)));
                    add_relation_between_new_and_ith_input(input, i, j, neuroLib.idtRange);
                }
                cout << " bir pushed" << endl;
            }
        }
    }
    input.assert_size();
}

unsigned int jbilinear_HH(Neuron &neuron, NeuroLib &neuroLib, Input &input, jND &jnd, Cross &cross, double end_t, double ignore_t, size corrSize, vector<double> &tsp, double pairs[], bool type[], double tau_er, double tau_ed, double tau_ir, double tau_id, double vC, double vB, int afterCrossBehavior, bool spikeShape) {
    size i, j, i_prior_cross;
    double tstep = neuroLib.tstep;
    double iend = end_t/tstep;
    double tol_tl = neuroLib.nt0;
    double tol_tb = neuroLib.idtRange[ndt-1] - round(ignore_t/tstep);
    std::cout << "linear corr length " << tol_tl << std::endl;
    std::cout << "bilinear corr length " << tol_tb << std::endl;
    std::cout << "total inputs " << neuron.tin.size() << std::endl;
    size tail_b = 0, tail_l = 0;
    size old_tail_b, old_tail_l;
    bool crossed, spiked;
    double tref = neuron.tref/tstep;
    size i_, j_;
    double r_;
    bool debug = false;
    size ii;
    jb::getNear(neuroLib.vRange, neuroLib.nv, cross.v[0], r_, i_, j_);
    cross.vCross.push_back(IJR(i_,j_,r_));
    jnd.t.push_back(0);
    jnd.v.push_back(cross.v[0]);
    for (i=0;i<neuron.tin.size();i++) {
        input.assert_size();
        input.t.push_back(neuron.tin[i]/tstep);
        assert(input.t.size() == i+1);
        old_tail_l = tail_l;
        old_tail_b = tail_b;
        move_corr_window(neuron.tin, tail_l, input.t[i], tol_tl, tstep);
        move_corr_window(neuron.tin, tail_b, input.t[i], tol_tb, tstep);
        crossed = update_vinit_of_new_input_check_crossing(input, cross, neuroLib, neuron, i, tail_l, tail_b, tol_tl, tol_tb, end_t, jnd, t_cross, vC, corrSize, debug);
        ii = 0;
        while (crossed) {
            ii++;
            if (t_cross < neuron.tin[i]/tstep + 1e-14) {
                cout << " unlikely, cross upon or before input" << endl;
                tail_l = old_tail_l;
                tail_b = old_tail_b;
                i_prior_cross = i-1;
            } else {
                i_prior_cross = i
            }
            if (jb::debug) {
                cout << " vCross = " << jnd.v.back();
                cout << " tCross = " << jnd.t.back() << endl;
                cout << " input " << i_prior_cross << ", t" << input.t[i_prior_cross]*tstep << endl;
            }
            if (input.t.size() != i_prior_cross+1) {
                cout << " input size not correct" << endl;
                assert(input.t.size() == i_prior_cross+1);
            }
            nc_old = nc;
            if (spikeShape) {
            } else {
                tBack = t_cross+tref;
                if (tBack > iend) {
                    tBack = iend;
                }
                vBack = neuron.vReset;
                double tmpTsp = t_cross*tstep + neuron.tref/2;
                if (tmpTsp <= end_t) {
                    tsp.push_back(tmpTsp);
                    nc++;
                }
                spiked = nc - nc_old;
                i = i_prior_cross + 1;
                if (i < neuron.inID.size()) {
                    while (neuron.tin[i]-1e-10 < tBack*tstep + 1e-14) {
                        i++;
                        if (i==neuron.inID.size()) {
                            break;
                        }
                    }
                }
                i--;
            }
            jnd.t.push_back(tBack);
            jnd.v.push_back(vBack);
            if (jb::debug) {
                cout << " backed at " << tBack*neuroLib.tstep << endl;
                cout << " input from " << i_prior_cross << " to " << i << endl;
                cout << " input at " << neuron.tin[i] << endl;
            }
            if (vBack >= vC) {
                break;
            }
            jb::getNear(neuroLib.vRange,neuroLib.nv,
                            vBack, r_, i_, j_);
            cross.vCross.push_back(IJR(i_,j_,r_));
            cross.nCross++;
            cross.iCross.push_back(cross.v.size());
            cross.tCross.push_back(tBack);
            cout << "crossed " << endl;

            old_tail_b = tail_b;
            old_tail_l = tail_l;
            move_corr_window(neuron.tin, tail_l, tBack, tol_tl, tstep);
            move_corr_window(neuron.tin, tail_b, tBack, tol_tb, tstep);
            update_info_after_cross(input, neuroLib, cross, neuron, tBack, vBack, i_prior_cross, tail_l, tail_b, i, corrSize, tsp, spiked, afterCrossBehavior);
            cout << " updated" << endl;
            if (!afterCrossBehavior) {
                tail_b = i + 1;
                tail_l = i + 1;
            } else {
                if (afterCrossBehavior == 1) {
                    tail_b = i + 1;
                }
            }
            crossed = check_crossing(input, neuroLib, cross, neuron, tol_tl, tol_tb, end_t, tail_l, tail_b, i, jnd, t_cross, debug);
            cout << "------------" << endl;
        }
        //cout << " dead with " << i << ", v " << jnd.v.back() << " at " << jnd.t.back() << endl;
        //input.print_this(i);
    }
    cout << " crossed " << cross.nCross << " times." << endl;
    return nc;
}
