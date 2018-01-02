function dh=gatingh(v,h,vt,t)
    if nargin <4
        dh=alphah(v,vt).*(1-h)-betah(v,vt).*h;
    else
        alpha = alphah(v,vt);
        beta = betah(v,vt);
        tau_h = 1./(alpha+beta);
        h_inf = alpha.*tau_h;
        tau = exp(-t./tau_h);
        dh = h_inf.*(1-tau) + h.*tau;
    end
end
