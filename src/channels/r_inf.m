function r_inf = r_inf(v)
    alpha = alphar(v);
    beta = betar(v);
    tau_r = 1./(alpha+beta);
    r_inf = alpha.*tau_r;
end