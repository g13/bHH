function y = tau_u(V,Vx)

    y=30.8+(211.4+exp((V+Vx+113.2)/5))./(3.7*(1+exp((V+Vx+84)/3.2)));
    %y
    %y=30.8+(211.4+exp((V+Vx)/5)*exp(113.2/5))./(3.7*(1+exp((V+Vx)/3.2)*exp(84/3.2)));
    %y

% may be a typo in the paper, 30.8 should be outside the fraction
% the original:
%  y=(30.8+211.4+exp((V+Vx+113.2)/5))./(3.7*(1+exp((V+Vx+84)/3.2)));
%  y
%  y=(30.8+211.4 + exp((V+Vx)/5) * exp(113.2/5) )/( 3.7*( 1 + exp((V+Vx)/3.2) * exp(84/3.2) ) );
%  y

end
