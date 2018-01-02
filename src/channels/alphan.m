function y = alphan(V,V_T)

    y=-0.032*(V-V_T-15)./(exp(-(V-V_T-15)/5)-1);
    idx = V==V_T+15;
    y(idx)=0.16; % inf in the denominator

end
