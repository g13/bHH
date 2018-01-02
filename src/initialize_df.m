function [v,m,n,h,p,q,r,s,u] = initialize_df(v,para,ipick)

v=v;
m=alpham(v,para.vT(ipick))./(alpham(v,para.vT(ipick))+betam(v,para.vT(ipick)));
n=alphan(v,para.vT(ipick))./(alphan(v,para.vT(ipick))+betan(v,para.vT(ipick)));
h=alphah(v,para.vT(ipick))./(alphah(v,para.vT(ipick))+betah(v,para.vT(ipick)));
p=p_inf(v);
q=alphaq(v)./(alphaq(v)+betaq(v));
r=alphar(v)./(alphar(v)+betar(v));
s=s_inf(v,para.vX(ipick));
u=u_inf(v,para.vX(ipick));

end
