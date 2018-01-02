function dv=volt_if(t,v,para,eif)
    dv = -para.gLeak.*(v-para.vLeak - eif*para.DeltaT.*exp((v-para.vT)./para.DeltaT));
    [ge,gi] = synapticCond(para,t);
    synapticCurrent = -sum([ge,gi].*[(v-para.vE),(v-para.vI)],2);
    dv = dv + synapticCurrent + para.current(t,para.fCurrent./para.S);
end
