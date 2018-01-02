#ifndef CHANNEL_H
#define CHANNEL_H
#include <cmath>
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

inline double fk_HH(double v, double m, double n, double h, double gE, double gI, double *pairs){
    // pairs: gNa,vNa,gK,vK,gLeak,vLeak,vE,vI
    return -pairs[0]*(v-pairs[1])*m*m*m*h - pairs[2]*n*n*n*n*(v-pairs[3]) - pairs[4]*(v-pairs[5]) - gE*(v-pairs[6]) - gI*(v-pairs[7]);
}

inline double fk_IF(double v, double gE, double gI, double *pairs, bool eif) {
    double dv = -pairs[4]*(v-pairs[5]) - gE*(v-pairs[6]) - gI*(v-pairs[7]);
    if (eif)
        return dv - pairs[4]*pairs[10]*exp((v-pairs[8])/pairs[10]);
    else
        return dv;
}
#endif
