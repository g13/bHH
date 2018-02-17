#ifndef CHANNEL_H
#define CHANNEL_H
#include <cmath>
inline double u_inf(double V, double Vx) {
    return 1.0/(1.0+exp((V+Vx+81.0)/4.0));
}
inline double tau_u(double V, double Vx) {
    return 30.8+(211.4+exp((V+Vx+113.2)/5.0))/(3.7*(1.0+exp((V+Vx+84.0)/3.2)));
}
inline double gatingu(double V, double u, double Vx) {
    return (u_inf(V,Vx)-u)/tau_u(V,Vx);
}
inline double s_inf(double V, double Vx) {
    return 1.0/(1.0+exp(-(V+Vx+57.0)/6.2));
}
inline double alphar(double V) {
    return 0.000457*exp((-13.0-V)/50.0);
}
inline double betar(double V) {
    return 0.0065/(exp((-15.0-V)/28.0)+1.0);
}
inline double gatingr(double V, double r) {
    return alphar(V)*(1.0-r) - betar(V)*r;
}
inline double alphaq(double V) {
    if (V == -27.0) {
        return 0.2090;
    } else {
        return 0.055*(-27.0-V)/(exp((-27.0-V)/3.8)-1.0);
    }
}
inline double betaq(double V) {
    return 0.94*exp((-75.0-V)/17.0);
}
inline double gatingq(double V, double q) {
    return alphaq(V)*(1.0-q)-betaq(V)*q;
}
inline double tau_p(double V, double tau_max){
    return tau_max/(3.3*exp((V+35.0)/20.0)+exp(-(V+35.0)/20.0));
}
inline double p_inf(double V) {
    return 1.0/(1.0+exp(-(V+35.0)/10.0));
}
inline double gatingp(double V, double p, double tau_max) {
    return (p_inf(V)-p)/tau_p(V,tau_max);
}
inline double alpham(double V, double V_T) {
    if (V-V_T != 13.0)
        return -0.32*(V-V_T-13.0)/(exp(-(V-V_T-13.0)/4.0)-1.0);
    else return 0.32*4.0;
}
inline double betam(double V, double V_T) {
    if (V-V_T != 40.0)
        return 0.28*(V-V_T-40.0)/(exp((V-V_T-40.0)/5.0)-1.0);
    else return 0.28*5.0;
}
inline double alphan(double V, double V_T) {
    if (V!=V_T+15.0)
        return -0.032*(V-V_T-15.0)/(exp(-(V-V_T-15.0)/5.0)-1.0);
    else return 0.16; 
}

inline double betan(double V, double V_T) {
    return 0.5*exp(-(V-V_T-10.0)/40.0);
}

inline double alphah(double V, double V_T) {
    return 0.128*exp(-(V-V_T-17.0)/18.0);
}
inline double betah(double V, double V_T) {
    return 4.0/(1.0+exp(-(V-V_T-40.0)/5.0));
}
inline double r_inf(double v){
    double alpha = alphar(v);
    double beta = betar(v);
    double tau_r = 1.0/(alpha+beta);
    return alpha*tau_r;
}
inline double q_inf(double v){
    double alpha = alphaq(v);
    double beta = betaq(v);
    double tau_q = 1.0/(alpha+beta);
    return alpha*tau_q;
}
inline double m_inf(double v, double vt){
    double alpha = alpham(v,vt);
    double beta = betam(v,vt);
    double tau_m = 1.0/(alpha+beta);
    return alpha*tau_m;
}
inline double n_inf(double v, double vt){
    double alpha = alphan(v,vt);
    double beta = betan(v,vt);
    double tau_n = 1.0/(alpha+beta);
    return alpha*tau_n;
}
inline double h_inf(double v, double vt){
    double alpha = alphah(v,vt);
    double beta = betah(v,vt);
    double tau_h = 1.0/(alpha+beta);
    return alpha*tau_h;
}
inline double gatingm(double v, double m, double vt) {
    double alpha = alpham(v,vt);
    return alpha-(betam(v,vt)+alpha)*m;
}
inline double gatingn(double v, double n, double vt) {
    double alpha = alphan(v,vt);
    return alpha-(betan(v,vt)+alpha)*n;
}
inline double gatingh(double v, double h, double vt) {
    return alphah(v,vt)*(1-h)-betah(v,vt)*h;
}

inline double fk_HH(double v, double m, double n, double h, double p, double q, double r, double s, double u, double gE, double gI, double *pairs, bool *type){
    // pairs: gNa,vNa,gK,vK,gLeak,vLeak,vE,vI,
    //         0   1   2  3  4    5      6  7 
    //        vT,vRest,DeltaT,
    //         8   9    10
    //        gM,gL,gT,vCA,Vx,tau_max
    //        11 12 13  14 15   16
    double currents = - pairs[4]*(v-pairs[5]) - gE*(v-pairs[6]) - gI*(v-pairs[7]);
    if (type[0])
        currents += -pairs[0]*(v-pairs[1])*m*m*m*h;
    if (type[1])
        currents += -pairs[2]*(v-pairs[3])*n*n*n*n;
    if (type[2])
        currents += -pairs[11]*(v-pairs[3])*p;
    if (type[3])
        currents += -pairs[12]*(v-pairs[14])*q*q*r;
    if (type[4])
        currents += -pairs[13]*(v-pairs[14])*s*s*u;
    return currents;
}

inline double fk_IF(double v, double gE, double gI, double *pairs, bool eif) {
    double dv = -pairs[4]*(v-pairs[5]) - gE*(v-pairs[6]) - gI*(v-pairs[7]);
    if (eif)
        return dv - pairs[4]*pairs[10]*exp((v-pairs[8])/pairs[10]);
    else
        return dv;
}
#endif
