function dv=volt_df(t,v,m,n,h,p,q,r,s,u,para,bool,ipick,nf)

    currents = [ -para.gNa(ipick).*m.^3.*h.*(v-para.vNa),...
                 -para.gK(ipick).*n.^4.*(v-para.vK),...
                 -para.gM(ipick).*p.*(v-para.vK),...
                 -para.gL(ipick).*q.^2.*r.*(v-para.vCA),...
                 -para.gT(ipick).*s.^2.*u.*(v-para.vCA)...
                ];
    idx = isnan(currents);
    currents(idx) = 0;
    idx = isinf(currents);
    currents(idx) = 0;
    [ge,gi] = synapticCond_df(para,t,ipick);
    synapticCurrent = -sum([ge,gi].*[(v-para.vE),(v-para.vI)],2);
    %dv = (-para.gLeak(ipick).*(v-para.vLeak(ipick)) + para.current(t,para.fCurrent./para.S(ipick)) + sum(currents.*repmat(bool(para.type(ipick),:),[nf,1]),2)) + synapticCurrent;
    dv = (-para.gLeak(ipick).*(v-para.vLeak(ipick)) + sum(currents.*repmat(bool(para.type(ipick),:),[nf,1]),2)) + synapticCurrent;
end
