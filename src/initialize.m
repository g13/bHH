function [v,m,n,h,p,q,r,s,u] = initialize(v,para)

v=v;
m=alpham(v,para.vT)./(alpham(v,para.vT)+betam(v,para.vT));
n=alphan(v,para.vT)./(alphan(v,para.vT)+betan(v,para.vT));
h=alphah(v,para.vT)./(alphah(v,para.vT)+betah(v,para.vT));
p=p_inf(v);
q=alphaq(v)./(alphaq(v)+betaq(v));
r=alphar(v)./(alphar(v)+betar(v));
s=s_inf(v,para.vX);
u=u_inf(v,para.vX);

end
