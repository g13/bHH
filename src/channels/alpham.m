function y = alpham(V,V_T) 

   y=-0.32*(V-V_T-13)./(exp(-(V-V_T-13)/4)-1);
   idx = V-V_T-13==0;
   y(idx) = 0.32*4;

end
