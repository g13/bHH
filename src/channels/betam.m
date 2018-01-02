function y = betam(V,V_T)

    y=0.28*(V-V_T-40)./(exp((V-V_T-40)/5)-1);
    idx = V-V_T == 40;
    y(idx) = 0.28*5;

end
