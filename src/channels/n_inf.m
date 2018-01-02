function n_inf = n_inf(v,vt)
    alpha = alphan(v,vt);
    beta = betan(v,vt);
    tau_n = 1./(alpha+beta);
    n_inf = alpha.*tau_n;
end