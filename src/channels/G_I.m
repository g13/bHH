function a = G_I(t0,t,f_I,Time_InConR,Time_InCon,gConI)

a=f_I*gConI*(exp(-(t-t0)/Time_InCon)-exp(-(t-t0)/Time_InConR));
idx = repmat((t-t0),size(f_I))<0;
a(idx) = 0;
end
