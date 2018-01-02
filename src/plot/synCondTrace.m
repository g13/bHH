function [gE,gI]= synCondTrace(para,t,j)
    if nargin < 3
        j = 1;
    end
    gE = zeros(size(t)); 
    for i=1:para.nEspike
        gE = gE + G_E(para.tE(i),t,para.f_E./para.S(j),para.tau_er,para.tau_e,para.gConE);
    end
    gI = zeros(size(t)); 
    for i=1:para.nIspike
        gI = gI + G_I(para.tI(i),t,para.f_I./para.S(j),para.tau_ir,para.tau_i,para.gConI);
    end
end