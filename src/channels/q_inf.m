function q_inf = q_inf(v)
    alpha = alphaq(v);
    beta = betaq(v);
    tau_q = 1./(alpha+beta);
    q_inf = alpha.*tau_q;
end