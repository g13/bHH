function y = alphaq(V)

    y=0.055*(-27-V)./(exp((-27-V)/3.8)-1);
    idx = -27-V==0;
    y(idx) = 0.2090;

end
    
    
