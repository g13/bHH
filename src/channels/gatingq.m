function dq=gatingq(v,q,t)
    if nargin < 3 
        dq=alphaq(v).*(1-q)-betaq(v).*q;
    else
        alpha = alphaq(v);
        beta = betaq(v);
        tau_q = 1./(alpha+beta);
        q_inf = alpha.*tau_q;
        tau = exp(-t./tau_q);
        dq = q_inf.*(1-tau) + q.*tau;
    end
end
