function a = G_E_df(t0,t,f_E,Time_ExConR,Time_ExCon,gConE)

a = f_E .* gConE.*(exp(-(t-t0)/Time_ExCon)-exp(-(t-t0)/Time_ExConR));
idx = (t-t0)<0;
a(idx) = 0;

end
