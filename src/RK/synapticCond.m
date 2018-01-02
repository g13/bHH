function [gE,gI]= synapticCond(para,t)
    gE = sum(G_E(para.tE,t,para.f_E./para.S,para.tau_er,para.tau_e,para.gConE),2);
    gI = sum(G_I(para.tI,t,para.f_I./para.S,para.tau_ir,para.tau_i,para.gConI),2);
end