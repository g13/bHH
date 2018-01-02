function dv=voltagev(t,v,m,n,h,p,q,r,s,u,para,bool)

    currents = [ -para.gNa.*m.^3.*h.*(v-para.vNa),...
                 -para.gK.*n.^4.*(v-para.vK),...
                 -para.gM.*p.*(v-para.vK),...
                 -para.gL.*q.^2.*r.*(v-para.vCA),...
                 -para.gT.*s.^2.*u.*(v-para.vCA)...
                ];
    idx = isnan(currents);
    currents(idx) = 0;
    idx = isinf(currents);
    currents(idx) = 0;
    dv = -para.gLeak.*(v-para.vLeak) + sum(currents.*bool(para.type,:),2);
    [ge,gi] = synapticCond(para,t);
    synapticCurrent = -sum([ge,gi].*[(v-para.vE),(v-para.vI)],2);
    dv = dv + synapticCurrent + para.current(t,para.fCurrent./para.S);
end
