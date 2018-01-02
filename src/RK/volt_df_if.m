function dv=volt_df_if(t,v,para,ipick,eif)
    dv = -para.gLeak(ipick).*(v-para.vLeak(ipick) - eif*para.DeltaT(ipick).*exp((v-para.vT(ipick))/para.DeltaT(ipick)));
    [ge,gi] = synapticCond_df(para,t,ipick);
    synapticCurrent = -sum([ge,gi].*[(v-para.vE),(v-para.vI)],2);
    dv = dv + synapticCurrent + para.current(t,para.fCurrent./para.S(ipick));
end
