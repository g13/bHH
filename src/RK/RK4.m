function y = RK4(name,v0,para,bool,N,tstep,etime,savedata)  %RK4 method
if nargin < 8
    savedata = false;
end
if isequal(name,'TR_somato_Rat');
    para.VK = -100;
end
stime=0;
count=round((etime-stime)./tstep)+1;

y=zeros(9,count,N);
t=stime;

[v,m,n,h,p,q,r,s,u] = initialize(v0,para);

yold = [v,m,n,h,p,q,r,s,u];
y(:,1,:)=yold';

K=zeros(N,9,4);
tstep2 = 0.5*tstep;
for i=1:count-1
    K(:,1,1)=voltagev(t,v,m,n,h,p,q,r,s,u,para,bool); v2 = v+tstep2*K(:,1,1);
    K(:,2,1)=gatingm(v,m,para.vT);                    m2 = m+tstep2*K(:,2,1); 
    K(:,3,1)=gatingn(v,n,para.vT);                    n2 = n+tstep2*K(:,3,1);
    K(:,4,1)=gatingh(v,h,para.vT);                    h2 = h+tstep2*K(:,4,1);
    K(:,5,1)=gatingp(v,p,para.tau_max);               p2 = p+tstep2*K(:,5,1);
    K(:,6,1)=gatingq(v,q);                            q2 = q+tstep2*K(:,6,1);
    K(:,7,1)=gatingr(v,r);                            r2 = r+tstep2*K(:,7,1);
    K(:,8,1)=  s_inf(v2,para.vX);                     s2 =          K(:,8,1);
    K(:,9,1)=gatingu(v,u,para.vX);                    u2 = u+tstep2*K(:,9,1);
    t2 = t+tstep2; 

    K(:,1,2)=voltagev(t2,v2,m2,n2,h2,p2,q2,r2,s2,u2,para,bool); v3 = v+tstep2*K(:,1,2);
    K(:,2,2)=gatingm(v2,m2,para.vT);                            m3 = m+tstep2*K(:,2,2); 
    K(:,3,2)=gatingn(v2,n2,para.vT);                            n3 = n+tstep2*K(:,3,2);
    K(:,4,2)=gatingh(v2,h2,para.vT);                            h3 = h+tstep2*K(:,4,2);
    K(:,5,2)=gatingp(v2,p2,para.tau_max);                       p3 = p+tstep2*K(:,5,2);
    K(:,6,2)=gatingq(v2,q2);                                    q3 = q+tstep2*K(:,6,2);
    K(:,7,2)=gatingr(v2,r2);                                    r3 = r+tstep2*K(:,7,2);
    K(:,8,2)=  s_inf(v3,para.vX);                               s3 =          K(:,8,2);
    K(:,9,2)=gatingu(v2,u2,para.vX);                            u3 = u+tstep2*K(:,9,2);

    K(:,1,3)=voltagev(t2,v3,m3,n3,h3,p3,q3,r3,s3,u3,para,bool); v4 = v+tstep*K(:,1,3);
    K(:,2,3)=gatingm(v3,m3,para.vT);                            m4 = m+tstep*K(:,2,3); 
    K(:,3,3)=gatingn(v3,n3,para.vT);                            n4 = n+tstep*K(:,3,3);
    K(:,4,3)=gatingh(v3,h3,para.vT);                            h4 = h+tstep*K(:,4,3);
    K(:,5,3)=gatingp(v3,p3,para.tau_max);                       p4 = p+tstep*K(:,5,3);
    K(:,6,3)=gatingq(v3,q3);                                    q4 = q+tstep*K(:,6,3);
    K(:,7,3)=gatingr(v3,r3);                                    r4 = r+tstep*K(:,7,3);
    K(:,8,3)=  s_inf(v4,para.vX);                               s4 =         K(:,8,3);
    K(:,9,3)=gatingu(v3,u3,para.vX);                            u4 = u+tstep*K(:,9,3);

    t = t+tstep;
    K(:,1,4)=voltagev(t,v4,m4,n4,h4,p4,q4,r4,s4,u4,para,bool);
    K(:,2,4)=gatingm(v4,m4,para.vT);
    K(:,3,4)=gatingn(v4,n4,para.vT);
    K(:,4,4)=gatingh(v4,h4,para.vT);
    K(:,5,4)=gatingp(v4,p4,para.tau_max);
    K(:,6,4)=gatingq(v4,q4);
    K(:,7,4)=gatingr(v4,r4);
    K(:,8,4)=  s_inf(v4,para.vX);
    K(:,9,4)=gatingu(v4,u4,para.vX);

    %y(:,i+1,:)=y(:,i,:)+tstep/6*(K(:,:,1)+2*K(:,:,2)+2*K(:,:,3)+K(:,:,4));
    ynew = yold + tstep/6*(K(:,:,1)+2*(K(:,:,2)+K(:,:,3))+K(:,:,4));
    yold = ynew;
    v=ynew(:,1);
    m=ynew(:,2);
    n=ynew(:,3);
    h=ynew(:,4);
    p=ynew(:,5);
    q=ynew(:,6);
    r=ynew(:,7);
    s=s_inf(v,para.vX);
    u=ynew(:,9);
    y(:,i+1,:) = [v,m,n,h,p,q,r,s,u]';
end 
if savedata
    save(['data-',name,'.mat'],'y');
end
end
