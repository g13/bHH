function dr=gatingr(v,r,t)
    if nargin < 3    
        dr=alphar(v).*(1-r)-betar(v).*r;
    else
        alpha = alphar(v);
        beta = betar(v);
        tau_r = 1./(alpha+beta);
        r_inf = alpha.*tau_r;
        tau = exp(-t./tau_r);
        dr = r_inf.*(1-tau) + r.*tau;
    end
end
 
