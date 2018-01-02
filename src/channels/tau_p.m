function y=tau_p(V,tau_max)

    y=tau_max./(3.3*exp((V+35)/20)+exp(-(V+35)/20));

end
