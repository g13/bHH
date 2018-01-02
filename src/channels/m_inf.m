function m_inf = m_inf(v,vt)
        alpha = alpham(v,vt);
        beta = betam(v,vt);
        tau_m = 1./(alpha+beta);
        m_inf = alpha.*tau_m;
end