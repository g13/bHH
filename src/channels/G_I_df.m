function a = G_I_df(t0,t,f_I,Time_InConR,Time_InCon,gConI)

a=f_I.*gConI.*(exp(-(t-t0)/Time_InCon)-exp(-(t-t0)/Time_InConR));
idx = (t-t0)<0;
a(idx) = 0;
end
