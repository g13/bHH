% scan input strengthes and record the corresponding EPSP, IPSP and SSP
% save the data for processing in polynomial_rule.m
function v0Iterate(tstep,dur,name,f_E,f_I,v0,vThresh,vorigin,theme)
    if nargin < 9
        theme = '';
        if nargin < 8
            vorigin = false;
            if nargin < 7
                vThresh = 0.0;
            end
        end
    end
    load(['parameters-',name],'para','bool','type','species','posInBrain','n');
    tcount=round(dur./tstep)+1;
    para.fCurrent = 0;
    nv0 = length(v0);
    v0 = para.vRest*ones(1,nv0) + ((para.vT - para.vRest)*v0);
    IPSP = zeros(tcount,nv0,n);
    EPSP = zeros(tcount,nv0,n);
    SSP = zeros(tcount,nv0,n);
    
    para.f_I = f_I;
    para.f_E = 0;
    for i=1:nv0
        x=RK4(name,v0(:,i),para,bool,n,tstep,dur);
        v = squeeze(x(1,:,:));
        IPSP(:,i,:)=v;
    end

    para.f_E = f_E;
    para.f_I = 0;
    for i=1:nv0
        x=RK4(name,v0(:,i),para,bool,n,tstep,dur);
        v = squeeze(x(1,:,:));
        EPSP(:,i,:)=v;
    end

    para.f_E = f_E;
    para.f_I = f_I;
    for i=1:nv0
        x=RK4(name,v0(:,i),para,bool,n,tstep,dur);
        v = squeeze(x(1,:,:));
        SSP(:,i,:)=v;
    end
    if vorigin
        vrest0 = reshape(para.vRest,[1,1,n]);
        vrest = repmat(vrest0,[tcount,nv0,1]);
    else
        vrest = zeros(tcount,nv0,n);
        para.f_E = 0.0;
        para.f_I = 0.0;
        for i=1:nv0
            x=RK4(name,v0(:,i),para,bool,n,tstep,dur);
            v = squeeze(x(1,:,:));
            vrest(:,i,:)=v;
        end
    end

    Emax = max(EPSP);
    gotEspike = Emax>vThresh; 
    gotEspike = reshape(gotEspike,[nv0,n]);
    Imax = max(IPSP);
    gotIspike = Imax>vThresh; 
    gotIspike = reshape(gotIspike,[nv0,n]);
    
    VI=IPSP-vrest;
    VI = reshape(VI,[tcount,nv0*n]);
    VI = mat2cell(VI,tcount, nv0*ones(n,1));

    VE=EPSP-vrest;
    VE = reshape(VE,[tcount,nv0*n]);
    VE = mat2cell(VE,tcount, nv0*ones(n,1));

    VS= SSP-vrest;
    VS = reshape(VS,[tcount,nv0*n]);
    VS = mat2cell(VS,tcount, nv0*ones(n,1));
    
    save(['data-v0-',name,theme,'.mat'],'tstep','VI','VE','VS','nv0','gotEspike','gotIspike');
end
