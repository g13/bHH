function y = RK2_IF(name,v0,para,bool,N,tstep,etime,savedata,eif,reset)  %RK4 method
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
    
    y=zeros(1,count,N);
    t=stime;
    
    v = v0;
    
    y(:,1,:)=v';
    if eif
        vT = para.vNa*ones(N,1);
    else
        vT = para.vT;
    end
    
    K=zeros(N,1,2);
    for i=1:count-1
        told = t;
        K(:,1,1)=volt_if(t,v,para,eif); v2 = v+tstep*K(:,1,1);
        t = t+tstep;
        K(:,1,2)=volt_if(t,v2,para,eif);
        vnew = v + tstep/2*(K(:,:,1)+K(:,:,2));
        if reset
            for j=1:N
                if vnew(j) > vT(j)
                    dtsp = (vT(j)-v(j))/(vnew(j)-v(j))*tstep;
                    K(:,1,1)=volt_if(told+dtsp,para.vLeak,para,eif); v2 = para.vLeak+(tstep-dtsp)*K(:,1,1);
                    K(:,1,2)=volt_if(t,v2,para,eif);
                    vnew(j) = para.vLeak(j) + (tstep-dtsp)/2*(K(j,1,1)+K(j,1,2));
                end
            end
        end
        v = vnew;
        y(1,i+1,:) = v';
    end 
    if savedata
        save(['data-',name,'.mat'],'y');
    end
end
