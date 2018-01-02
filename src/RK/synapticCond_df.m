function [gE,gI]= synapticCond_df(para,t,ipick)
    gE = sum(G_E_df(para.tE,t,para.f_E./para.S(ipick),para.tau_er,para.tau_e,para.gConE),2);
    gI = sum(G_I_df(para.tI,t,para.f_I./para.S(ipick),para.tau_ir,para.tau_i,para.gConI),2);
end
