function h_inf = h_inf(v,vt)
    alpha = alphah(v,vt);
    beta = betah(v,vt);
    tau_h = 1./(alpha+beta);
    h_inf = alpha.*tau_h;
end