function y = RK4_df(name,v0,para,bool,tstep,etime,ipick,savedata,gv0)  %RK4 method
if nargin < 9
    gv0.exist = false;
    if nargin < 8
        savedata = false;
    end
end
stime=0;
count=round((etime-stime)./tstep)+1;
% if isequal(size(para.f_E),size(para.f_I)) || (size(para.f_I,1) == 1 && sum(para.f_I)~=0) || (size(para.f_E,1) == 1 && sum(para.f_E)~=0)
%     nf = max(size(para.f_E,1),size(para.f_I,1));
%     if size(para.f_I,1) < size(para.f_E,1)
%         para.f_I = ones(nf,1)*para.f_I;
%     end
%     if size(para.f_E,1) < size(para.f_I,1)
%         para.f_E = ones(nf,1)*para.f_E;
%     end
%     if ~isequal(size(para.tI),size(para.f_I))
%         para.tI = ones(nf,1)*para.tI;
%     end
%     if ~isequal(size(para.tE),size(para.f_E))
%         para.tE = ones(nf,1)*para.tE;
%     end
% else
    assert(isequal(size(para.f_E),size(para.tE)) && isequal(size(para.f_I),size(para.tI)));
    nf = max(size(para.f_E,1),size(para.f_I,1));
% end


y=zeros(9,count,nf);
t=stime;

if para.vtime ~= -1
    v0 = ones(nf,1)*para.newv;
end

if ~gv0.exist
    [v,m,n,h,p,q,r,s,u] = initialize_df(v0,para,ipick);
else
    v = v0;
    if size(gv0.m0)~=size(v0)
        m = gv0.m0*ones(size(v0));
    else
        m = gv0.m0;
    end
    if size(gv0.n0)~=size(v0)
        n = gv0.n0*ones(size(v0));
    else
        n = gv0.n0;
    end
    if size(gv0.h0)~=size(v0)
        h = gv0.h0*ones(size(v0));
    else
        h = gv0.h0;
    end
    if size(gv0.p0)~=size(v0)
        p = gv0.p0*ones(size(v0));
    else
        p = gv0.p0;
    end
    if size(gv0.q0)~=size(v0)
        q = gv0.q0*ones(size(v0));
    else
        q = gv0.q0;
    end
    if size(gv0.r0)~=size(v0)
        r = gv0.r0*ones(size(v0));
    else
        r = gv0.r0;
    end
    if size(gv0.s0)~=size(v0)
        s = gv0.s0*ones(size(v0));
    else
        s = gv0.s0;
    end
    if size(gv0.u0)~=size(v0)
        u = gv0.u0*ones(size(v0));
    else
        u = gv0.u0;
    end
end

yold = [v,m,n,h,p,q,r,s,u];
y(:,1,:)=yold';

K=zeros(nf,9,4);
tstep2 = 0.5*tstep;
for i=1:count-1
    if round(para.vtime/tstep) + 1 < i+1
        K(:,1,1)=volt_df(t,v,m,n,h,p,q,r,s,u,para,bool,ipick,nf); v2 = v+tstep2*K(:,1,1);
        K(:,2,1)=gatingm(v,m,para.vT(ipick));                     m2 = m+tstep2*K(:,2,1); 
        K(:,3,1)=gatingn(v,n,para.vT(ipick));                     n2 = n+tstep2*K(:,3,1);
        K(:,4,1)=gatingh(v,h,para.vT(ipick));                     h2 = h+tstep2*K(:,4,1);
        K(:,5,1)=gatingp(v,p,para.tau_max(ipick));                p2 = p+tstep2*K(:,5,1);
        K(:,6,1)=gatingq(v,q);                                    q2 = q+tstep2*K(:,6,1);
        K(:,7,1)=gatingr(v,r);                                    r2 = r+tstep2*K(:,7,1);
        K(:,8,1)=  s_inf(v2,para.vX(ipick));                      s2 =          K(:,8,1);
        K(:,9,1)=gatingu(v,u,para.vX(ipick));                     u2 = u+tstep2*K(:,9,1);
        t2 = t+tstep2; 

        K(:,1,2)=volt_df(t2,v2,m2,n2,h2,p2,q2,r2,s2,u2,para,bool,ipick,nf); v3 = v+tstep2*K(:,1,2);
        K(:,2,2)=gatingm(v2,m2,para.vT(ipick));                             m3 = m+tstep2*K(:,2,2); 
        K(:,3,2)=gatingn(v2,n2,para.vT(ipick));                             n3 = n+tstep2*K(:,3,2);
        K(:,4,2)=gatingh(v2,h2,para.vT(ipick));                             h3 = h+tstep2*K(:,4,2);
        K(:,5,2)=gatingp(v2,p2,para.tau_max(ipick));                        p3 = p+tstep2*K(:,5,2);
        K(:,6,2)=gatingq(v2,q2);                                            q3 = q+tstep2*K(:,6,2);
        K(:,7,2)=gatingr(v2,r2);                                            r3 = r+tstep2*K(:,7,2);
        K(:,8,2)=  s_inf(v3,para.vX(ipick));                                s3 =          K(:,8,2);
        K(:,9,2)=gatingu(v2,u2,para.vX(ipick));                             u3 = u+tstep2*K(:,9,2);

        K(:,1,3)=volt_df(t2,v3,m3,n3,h3,p3,q3,r3,s3,u3,para,bool,ipick,nf); v4 = v+tstep*K(:,1,3);
        K(:,2,3)=gatingm(v3,m3,para.vT(ipick));                             m4 = m+tstep*K(:,2,3); 
        K(:,3,3)=gatingn(v3,n3,para.vT(ipick));                             n4 = n+tstep*K(:,3,3);
        K(:,4,3)=gatingh(v3,h3,para.vT(ipick));                             h4 = h+tstep*K(:,4,3);
        K(:,5,3)=gatingp(v3,p3,para.tau_max(ipick));                        p4 = p+tstep*K(:,5,3);
        K(:,6,3)=gatingq(v3,q3);                                            q4 = q+tstep*K(:,6,3);
        K(:,7,3)=gatingr(v3,r3);                                            r4 = r+tstep*K(:,7,3);
        K(:,8,3)=  s_inf(v4,para.vX(ipick));                                s4 =         K(:,8,3);
        K(:,9,3)=gatingu(v3,u3,para.vX(ipick));                             u4 = u+tstep*K(:,9,3);

        t = t+tstep;
        K(:,1,4)=volt_df(t,v4,m4,n4,h4,p4,q4,r4,s4,u4,para,bool,ipick,nf);
        K(:,2,4)=gatingm(v4,m4,para.vT(ipick));
        K(:,3,4)=gatingn(v4,n4,para.vT(ipick));
        K(:,4,4)=gatingh(v4,h4,para.vT(ipick));
        K(:,5,4)=gatingp(v4,p4,para.tau_max(ipick));
        K(:,6,4)=gatingq(v4,q4);
        K(:,7,4)=gatingr(v4,r4);
        K(:,8,4)=  s_inf(v4,para.vX(ipick));
        K(:,9,4)=gatingu(v4,u4,para.vX(ipick));

        %y(:,i+1,:)=y(:,i,:)+tstep/6*(K(:,:,1)+2*K(:,:,2)+2*K(:,:,3)+K(:,:,4));
        ynew = yold + tstep/6*(K(:,:,1)+2*(K(:,:,2)+K(:,:,3))+K(:,:,4));
        yold = ynew;
        v = ynew(:,1);
        m=ynew(:,2);
        n=ynew(:,3);
        h=ynew(:,4);
        p=ynew(:,5);
        q=ynew(:,6);
        r=ynew(:,7);
        s=s_inf(v,para.vX(ipick));
        u=ynew(:,9);
    else
        t = t+tstep;
    end
    y(:,i+1,:) = [v,m,n,h,p,q,r,s,u]';
end 
if savedata
    save(['data-',name,'.mat'],'y');
end
end
