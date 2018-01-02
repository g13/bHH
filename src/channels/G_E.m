function a = G_E(t0,t,f_E,Time_ExConR,Time_ExCon,gConE)

a = f_E * gConE*(exp(-(t-t0)/Time_ExCon)-exp(-(t-t0)/Time_ExConR));
idx = repmat((t-t0),size(f_E))<0;
a(idx) = 0;

end
