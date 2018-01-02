function [y,spiked] = RK2_IF_df(name,v0,para,bool,tstep,etime,ipick,savedata,eif,reset)  %RK4 method
    if nargin < 10
        reset = true;
        if nargin < 9
            eif = false;
            if nargin < 8
                savedata = false;
            end
        end
    end
    stime=0;
    count=round((etime-stime)./tstep)+1;
    assert(isequal(size(para.f_E),size(para.tE)) && isequal(size(para.f_I),size(para.tI)));
    nf = max(size(para.f_E,1),size(para.f_I,1));
    
    
    y=zeros(1,count,nf);
    t=stime;
    
    
    if round(para.vtime/tstep)+1 == 1
        v0 = ones(nf,1)*para.newv;
    end
    v = v0;
    y(:,1,:)=v';
    if eif
        vT = para.vNa(ipick);
    else
        vT = para.vT(ipick);
    end
    vR = para.vLeak(ipick)*ones(nf,1);
    K=zeros(nf,1,2);
    spiked = false;
    for i=1:count-1
        told = t;
        K(:,1,1)=volt_df_if(t,v,para,ipick,eif); v2 = v+tstep*K(:,1,1);
        t = t+tstep;
        K(:,1,2)=volt_df_if(t,v2,para,ipick,eif);
        vnew = v + tstep/2*(K(:,:,1)+K(:,:,2));
        for j=1:nf
            if vnew(j) > vT && reset
                dtsp = (vT-v(j))/(vnew(j)-v(j))*tstep;
                K(:,1,1)=volt_df_if(told+dtsp,vR,para,ipick,eif); v2 = vR+(tstep-dtsp)*K(:,1,1);
                K(:,1,2)=volt_df_if(t,v2,para,ipick,eif);
                vnew(j) = vR(j) + (tstep-dtsp)/2*(K(j,1,1)+K(j,1,2));
                spiked = true;
            end
        end
        if round(para.vtime/tstep) + 1 == i+1
            v = ones(nf,1)*para.newv;
        else
            v = vnew;
        end
        y(1,i+1,:) = v';
    end 
    if savedata
        save(['data-',name,'.mat'],'y');
    end
end
