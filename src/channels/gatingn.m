function dn=gatingn(v,n,vt,t)
    if nargin < 4 
        dn=alphan(v,vt).*(1-n)-betan(v,vt).*n;
    else
        alpha = alphan(v,vt);
        beta = betan(v,vt);
        tau_n = 1./(alpha+beta);
        n_inf = alpha.*tau_n;
        tau = exp(-t./tau_n);
        dn = n_inf.*(1-tau) + n.*tau;
    end
end
