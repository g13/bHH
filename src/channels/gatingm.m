function dm=gatingm(v,m,vt,t)
    if nargin < 4
        dm=alpham(v,vt).*(1-m)-betam(v,vt).*m;
    else
        alpha = alpham(v,vt);
        beta = betam(v,vt);
        tau_m = 1./(alpha+beta);
        m_inf = alpha.*tau_m;
        tau = exp(-t./tau_m);
        dm = m_inf.*(1-tau) + m.*tau;
    end
end
