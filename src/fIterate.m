% scan input strengthes and record the corresponding EPSP, IPSP and SSP
% save the data for processing in polynomial_rule.m
function fIterate(tstep,dur,name,f_E,f_I,vThresh,theme)
    if nargin < 7
        theme = '';
        if nargin < 6
            vThresh = 0.0;
        end
    end
    load(['parameters-',name],'para','bool','type','species','posInBrain','n');
    tcount=round(dur./tstep)+1;
    para.fCurrent = 0;
    v0 = para.vRest;
%     v0 = para.vLeak;
%     v0 = (para.vLeak + para.vT)/2;
    nI = length(f_I);
%     f_I = linspace(5e-6,5e-5,nI);
    IPSP = zeros(tcount,nI,n);
    nE = length(f_E);
%     f_E = linspace(3e-7,1e-5,nE);
    EPSP = zeros(tcount,nE,n);
    SSP = zeros(tcount,nI*nE,n);
    
    para.f_E = 0;
    for i=1:nI
        para.f_I = f_I(i);
        x=RK4(name,v0,para,bool,n,tstep,dur);
        v = squeeze(x(1,:,:));
        IPSP(:,i,:)=v;
    end
    
    para.f_I = 0;
    for i=1:nE
        para.f_E = f_E(i);
        x=RK4(name,v0,para,bool,n,tstep,dur);
        v = squeeze(x(1,:,:));
        EPSP(:,i,:)=v;
    end
     
    for i=1:nE
        para.f_E = f_E(i);
        for j=1:nI
            para.f_I = f_I(j);
            ind = (i-1)*nI+j;
            x=RK4(name,v0,para,bool,n,tstep,dur);
            v = squeeze(x(1,:,:));
            SSP(:,ind,:)=v;
         end
    end
    Emax = max(EPSP);
    gotEspike = Emax>vThresh; 
    gotEspike = reshape(gotEspike,[nE,n]);
    Imax = max(IPSP);
    gotIspike = Imax>vThresh; 
    gotIspike = reshape(gotIspike,[nI,n]);
    
    vrest0 = reshape(para.vRest,[1,1,n]);
    vrest = repmat(vrest0,[tcount,nI,1]);
    VI=IPSP-vrest;
    VI = reshape(VI,[tcount,nI*n]);
    VI = mat2cell(VI,tcount, nI*ones(n,1));

    vrest = repmat(vrest0,[tcount,nE,1]);
    VE=EPSP-vrest;
    VE = reshape(VE,[tcount,nE*n]);
    VE = mat2cell(VE,tcount, nE*ones(n,1));

    vrest = repmat(vrest0,[tcount,nI*nE,1]);
    VS= SSP-vrest;
    VS = reshape(VS,[tcount,nE*nI*n]);
    VS = mat2cell(VS,tcount, nE*nI*ones(n,1));
    
    save(['data-',name,theme,'.mat'],'tstep','VI','VE','VS','nE','nI','gotEspike','gotIspike');
end
