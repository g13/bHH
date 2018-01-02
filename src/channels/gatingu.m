function du = gatingu(v,u,Vx,t)
    if nargin < 4
        tau = tau_u(v,Vx);
        du=(u_inf(v,Vx)-u)./tau;
    else
        tau = exp(-t./tau_u(v,Vx));
        du = u_inf(v,Vx).*(1-tau) + u.*tau;
    end
end
