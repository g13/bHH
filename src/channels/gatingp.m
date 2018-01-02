function dp=gatingp(v,p,tau_max,t)
    if nargin < 4
        tau = tau_p(v,tau_max);
        dp=(p_inf(v)-p)./tau;
    else
        tau = exp(-t./tau_p(v,tau_max));
        dp = p_inf(v).*(1-tau) + p.*tau;
    end
end
