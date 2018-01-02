function [sEPSP,sIPSP,t] = noAdapV_k4(theme,name,pick,model,picformat,draw,ppp,loadData,npool,v0)
    global nvplot iv0Case;
    msgID = 'MATLAB:rankDeficientMatrix';
    warning('off',msgID);
    visible = true;
    if ~visible
        set(0,'DefaultFigureVisible','off');
    end
    FontSize = 16;
    set(0,'DefaultAxesFontSize',FontSize);
    set(0,'DefaultTextFontSize',FontSize-2);
    tstep = 1/32;
    l0 = 220; % estimating length
    dur = 300; % bilinear k length
    dur0 = 300; % linear length
    dtRange = [0:2:12,15:5:30,50,70,110,150,190,270];
    %dtRange = [0:6:12,18:12:30,54,78,126,174,212];
    %dtRange = 0:2:l0;
    ndt = length(dtRange);
    idtCase = [ndt-3,round(ndt/2),3];
    if nargin < 9
        v0 = -0.4:0.1:1.2;
        if nargin < 8
            npool = 1;
            if nargin < 7
                loadData = true;
                if nargin < 6
                    ppp = false;
                    if nargin < 5
                        draw = false;
                        if nargin <4
                            picformat = '';
                            if nargin < 3
                                model = 'HH';
                                if nargin < 2
                                    pick = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    %  paper
    width = 17;     height = width/16*9;
    % figure
    mleft = 0.0;    left = mleft*width;
    mbot = 0.0;     bot = mbot*height;    
    mright = 0.0;   fwidth = width-mleft-width*mright;
    mtop = 0.0;     fheight = height-mbot-height*mtop;
    pos0 = [width, height, left,bot,fwidth,fheight];  
    if ~isempty(picformat)
        if strcmp(picformat,'psc2')
            printDriver = ['-de',picformat];
            picformat = 'eps';
        else
            printDriver = ['-d',picformat];
        end
        dpi = '-r300';
    else
        dpi = '';
        printDriver = '';
    end
    switch model
        case 'HH'
            silico = @RK4_df;
        case 'IF'
%             silico = @RK2_IF_df;
            silico = @(a1,b2,c3,d4,e5,f6,g7,h8) RK2_IF_df(a1,b2,c3,d4,e5,f6,g7,h8,false,false);
        case 'EIF'
            silico = @RK2_EIF_df;
        otherwise
            disp('model not implemented');
            return
    end
    userpath([pwd,'/channels']);
    pname = ['parameters-',name,'-',theme,'-noAdap','-',model];
    disp(name);
    load([pname]);
    i = pick;
    ei = para.ei(i);
    dir = [name,'-',theme,'-',model];
    name = dir;
    if ~exist(dir,'dir')
        mkdir(dir);
    end
%     loadData = false;
    testEE = true;
    testII = true;
    testEI = true;
    testIE = true;
    singleStored = false;
    %singleStored = true;
%    multipleInput = false;
    multipleInput = true;
%    simpleTest = true;
    simpleTest = false;
    test = false;
%    test = true;
    assert(dtRange(ndt)<dur);
    assert((l0-dtRange(ndt-1))<=(dur-dtRange(ndt)));
    %idtCase = [1,round(ndt/2),ndt-1];
    ndtplot = length(idtCase);

    seed = 89749+1;
    ignore = 80;
    ignorefE = 0;
    ignorefI = 0; 
    rateE = 100/1000; % Hz/1000
    rateI = 120/1000;
    para.fCurrent = 0;
    %fE = linspace(0.0,1.0,5) * 1e-5;
    %fE = fE(2:5);
    %fE = (0.025:0.025:0.125) * 1e-6;
    %fE = (0.25:0.25:1.0) * 1e-5;
    %fE = (0.5:0.5:2.0) * 1e-5;
    fE = (0.125:0.125:0.5) * 1e-5;
    fE = fE(1:4);
    fI = (0.5:0.5:2.0) * 1e-5;
    fI = fI(1:4);
    %fI = fE;
%     fE = (0.5:0.5:1.0) * 1e-6;
%     fI = 0.5*1e-6;
    v0id = find(abs(v0 - 0.0)<1e-14);
    assert(~isempty(v0id));
    nv0 = length(v0);
    iv0Case = [1,round(nv0/2),nv0];
    nvplot = length(iv0Case);
    if ~loadData
        if npool > 1
            poolobj = gcp('nocreate'); % If no pool, do not create new one.
            if ~isempty(poolobj)
                if poolobj.NumWorkers ~= npool 
                    delete(poolobj);
                    parpool(npool);
                end
            else
                parpool(npool);
            end
        end
    else
        delete(gcp('nocreate'));
    end
    assert(nv0>1);
    nE = length(fE);
    nI = length(fI);
    tp0 = round(22/tstep);
    t = 0:tstep:dur;
    t0 = 0:tstep:dur0;
    nt = round(dur/tstep)+1;
    nt0 = round(dur0/tstep)+1;
    % 1 2 3
    % E I EI
    diri = [dir,'/',num2str(i)];
    if ~exist(diri,'dir')
        mkdir(diri);
    end
    vRange = para.vRest(i) + (para.vT(i)-para.vRest(i))*v0;
    if ~loadData
        if ~singleStored
            tic;
            [sEPSP,sIPSP,E_tmax,I_tmax,vleakage] = sPSP_check(silico,tstep,vRange,fI',fE',para,bool,name,dur,i,v0id,dtRange);
            toc;
            tic;
            if dur==dur0
                sEPSP0 = squeeze(sEPSP(:,:,1,:));
                sIPSP0 = squeeze(sIPSP(:,:,1,:));
                vleakage0 = squeeze(vleakage(:,1,:));
                h = plotsPSP0(sEPSP0,sIPSP0,vleakage0,fE,fI,nv0,dur,nt,v0id);
            else
                [h,sEPSP0,sIPSP0,vleakage0] = sPSP0_check(silico,tstep,vRange,fI',fE',para,bool,name,dur0,i,v0id);
            end
            delete([diri,'/*']);
            if draw
                fname = ['sPSP','-',name,'-',num2str(i)];
                printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
            end
            save(['single-',name,'-',num2str(i),'th.mat'],'sEPSP','sIPSP','vleakage','sEPSP0','sIPSP0','E_tmax','I_tmax','vleakage0');
            toc;
        else
            load(['single-',name,'-',num2str(i),'th.mat']);
            delete([diri,'/*']);
            if dur==dur0
                h = plotsPSP0(sEPSP0,sIPSP0,vleakage0,fE,fI,nv0,dur,nt,v0id);
            end
        end
        
        if ~strcmp(model,'HH')
            inv0 = nv0;
            para.vtime = -1;
            [~,tEmax] = max(sEPSP(:,nE,inv0));
            %m = m+vleakage(tEmax,inv0);
            para.f_I=0;
            para.tI=0;
            para.f_E=[fE(nE),fE(nE)];
            para.tE=[0,tEmax*tstep];
            [tmp,m1] = silico(name,vRange(inv0),para,bool,tstep,dur,i,false);
            para.tE=[0,0];
            [tmp,m2] = silico(name,vRange(inv0),para,bool,tstep,dur,i,false);
            while (m1 || m2) && inv0 > 1
            %while m>para.vRest(i)+v0(nv0)*2*(para.vT(i)-para.vRest(i)) && inv0 > 1
                disp(['spiked at ',num2str(inv0),'th v0']);
                inv0 = inv0-1;
                [~,tEmax] = max(sEPSP(:,nE,inv0));
                %m = m+vleakage(tEmax,inv0);
                para.tE=[0,tEmax*tstep];
                [tmp,m1] = silico(name,vRange(inv0),para,bool,tstep,dur,i,false);
                para.tE=[0,0];
                [tmp,m2] = silico(name,vRange(inv0),para,bool,tstep,dur,i,false);
            end
            disp(['subthreshold from ',num2str(inv0),'th v0'])
            disp(['total ',num2str(nv0),' v0']);
            assert(nv0 == inv0);
            %return
        end
        copyfile('noAdapV.m',[dir,'/noAdapV-',name,'-',num2str(i),'th.m']);
        disp(datestr(datetime('now')));
        %
        if draw
            h = figure;
            nsEPSP = sEPSP0./repmat(max(abs(sEPSP0)),[nt0,1]);
            nsIPSP = sIPSP0./repmat(max(abs(sIPSP0)),[nt0,1]);
            subplot(1,2,1)
            hold on
            c1 = linspace(0,5/6,nv0);
            c2 = linspace(0.3,1,nE)';
            c3 = 0.8;
            for iF = 1:nE
                for iv0 = 1:nv0
                    plot(t0,squeeze(nsEPSP(:,iF,iv0)),'Color',hsv2rgb([c1(iv0),c2(iF),c3]));
                end
            end
            title('norm. EPSP, hue v0, satur f');
            xlim([0,dur]);
            subplot(1,2,2)
            hold on
            c1 = linspace(0,5/6,nv0);
            c2 = linspace(0.3,1,nI)';
            c3 = 0.8;
            for iF = 1:nI
                for iv0 = 1:nv0
                    plot(t0,squeeze(nsIPSP(:,iF,iv0)),'Color',hsv2rgb([c1(iv0),c2(iF),c3]));
                end
            end
            title('norm. IPSP, hue v0, satur f');
            xlim([0,dur]);
            fname = ['nsPSP','-',name,'-',num2str(i)];
            printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
        end
        %
        
        % taking k at max |response|
%         [~,tQE] = max(sEPSP);
%         [~,tQI] = max(sIPSP);
%         tQE = round(mean(squeeze(tQE)));
%         tQI = round(mean(squeeze(tQI)));

        kVEE  = zeros(nt,ndt,ndt,nv0);
        kVII  = zeros(nt,ndt,ndt,nv0);
        kVEI  = zeros(nt,ndt,ndt,nv0);
        kVIE  = zeros(nt,ndt,ndt,nv0);

        if draw
            if testEE
                hEE0 = gobjects(1,ndtplot);
                for idt = 1:ndtplot
                    hEE0(idt) = figure;
                end
            end
            if testII
                hII0 = gobjects(1,ndtplot);
                for idt = 1:ndtplot
                    hII0(idt) = figure;
                end
            end
            if testEI
                hEI0 = gobjects(1,ndtplot);
                for idt = 1:ndtplot
                    hEI0(idt) = figure;
                end
            end
            if testIE
                hIE0 = gobjects(1,ndtplot);
                for idt = 1:ndtplot
                    hIE0(idt) = figure;
                end
            end
        end
        % redundant definition ignore
%         pv0 = zeros(nE*nI,nt);
%         vadd0 = zeros(nE*nI,nt);
        % % %
        idtplot = 0;
        incpercent = 10;
        ipercent = 1;
        tic;
        for idt = 1:ndt
            iidt = round(dtRange(idt)/tstep)+1;
            percent = idt/ndt*100;
            if floor(percent/incpercent) > ipercent
                disp(['--',num2str(round(percent)),'%dt']);
                disp(datestr(datetime('now')));
                ipercent = floor(percent/incpercent);
            end
            if sum((idt-idtCase)==0)==1 && draw
                dtplot = true;
                idtplot = idtplot+1;
            else
                dtplot = false;
            end
            
            if dtplot && ppp
                pp0 = true;
            else
                pp0 = false;
            end
            if testEE
                % EE
                vv0 = vRange(v0id)*ones(nE*nE,1);
                para.tE = repmat([0,dtRange(idt)],[nE*nE,1]);
                para.f_E = [reshape(ones(nE,1)*fE,[nE*nE,1]),repmat(fE',[nE,1])];
                para.tI = zeros(size(para.tE));
                para.f_I = zeros(size(para.f_E));
                [kVEE(:,:,idt,:),pv,vadd,vEEDoublet,vEE,h] = doubleCheck(silico,name,para,vv0,bool,tstep,iidt,dur,i,vleakage,sEPSP,sEPSP0,vRange,pp0,idt,ndt,dtRange);
                if dtplot
                    hEE0(idtplot) = drawExample(hEE0(idtplot),pv,vadd,vEEDoublet,squeeze(vleakage(:,idt,:)),t,iidt,tp0,kVEE(:,1,idt,:),nE*nE,vRange,para.vRest(i),vEE,'V_{E}V_{E}','EE');
                end
                if pp0
                    fname = [num2str(idt),'dt-EE'];
                    printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
                end
            end
            if testII
                % II
                vv0 = vRange(v0id)*ones(nI*nI,1);
                para.tI = repmat([0,dtRange(idt)],[nI*nI,1]);
                para.f_I = [reshape(ones(nI,1)*fI,[nI*nI,1]),repmat(fI',[nI,1])];
                para.tE = zeros(size(para.tI));
                para.f_E = zeros(size(para.f_I));
                [kVII(:,:,idt,:),pv,vadd,vIIDoublet,vII,h] = doubleCheck(silico,name,para,vv0,bool,tstep,iidt,dur,i,vleakage,sIPSP,sIPSP0,vRange,pp0,idt,ndt,dtRange);
                if dtplot
                    hII0(idtplot) = drawExample(hII0(idtplot),pv,vadd,vIIDoublet,squeeze(vleakage(:,idt,:)),t,iidt,tp0,kVII(:,1,idt,:),nI*nI,vRange,para.vRest(i),vII,'V_{I}V_{I}','II');
                end
                if pp0
                    fname = [num2str(idt),'dt-II'];
                    printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
                end
            end
            if testEI
                % EI 
                vv0 = vRange(v0id)*ones(nE*nI,1);
                para.f_E = reshape(ones(nI,1)*fE,[nE*nI,1]);
                para.tE = zeros(nE*nI,1);
                para.f_I = repmat(fI',[nE,1]);
                para.tI = dtRange(idt)*ones(nE*nI,1);
                [kVEI(:,:,idt,:),pv,vadd,vEIDoublet,vEI,h] = doubleCheck(silico,name,para,vv0,bool,tstep,iidt,dur,i,vleakage,sEPSP,sIPSP0,vRange,pp0,idt,ndt,dtRange);
                if dtplot
                    hEI0(idtplot) = drawExample(hEI0(idtplot),pv,vadd,vEIDoublet,squeeze(vleakage(:,idt,:)),t,iidt,tp0,kVEI(:,1,idt,:),nE*nI,vRange,para.vRest(i),vEI,'V_{E}V_{I}','EI');
                end
                if pp0
                    fname = [num2str(idt),'dt-EI'];
                    printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
                end
            end
            if testIE
                % IE
                vv0 = vRange(v0id)*ones(nE*nI,1);
                para.f_I = reshape(ones(nE,1)*fI,[nE*nI,1]);
                para.tI = zeros(nE*nI,1);
                para.f_E = repmat(fE',[nI,1]);
                para.tE = dtRange(idt)*ones(nE*nI,1);
                [kVIE(:,:,idt,:),pv,vadd,vIEDoublet,vIE,h] = doubleCheck(silico,name,para,vv0,bool,tstep,iidt,dur,i,vleakage,sIPSP,sEPSP0,vRange,pp0,idt,ndt,dtRange);
                if dtplot
                    hIE0(idtplot) = drawExample(hIE0(idtplot),pv,vadd,vIEDoublet,squeeze(vleakage(:,idt,:)),t,iidt,tp0,kVIE(:,1,idt,:),nE*nI,vRange,para.vRest(i),vIE,'V_{I}V_{E}','IE');
                end
                if pp0
                    fname = [num2str(idt),'dt-IE'];
                    printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
                end
            end
        end
        toc;
        disp('kV generated');
        save([name,'-',num2str(i),'th'],'kVEE','kVII','kVEI','kVIE','tp0','vleakage','sEPSP','sIPSP','dur','vRange','fE','fI','nE','nI','i','dtRange','sEPSP0','sIPSP0','vleakage0','ei','l0','tstep','dir','nt0','E_tmax','I_tmax');
        if draw
            assert(idtplot==ndtplot);
        end
    else
        if multipleInput
            load([name,'-',num2str(i),'th']);
            disp('kV loaded');
        end
    end
    if multipleInput
        kVEE = squeeze(kVEE(:,1,:,:));
        kVII = squeeze(kVII(:,1,:,:));
        kVEI = squeeze(kVEI(:,1,:,:));
        kVIE = squeeze(kVIE(:,1,:,:));
        disp('testing multiple inputs');
        % uniform K
            %>2 input
            % poisson inputs
        %  seed = 237587; % spike
%         poinumber = true;
        poinumber = false;
        poiend = false;
        idtRange = round(dtRange/tstep)+1;
        durpsp = l0-ignore;
%         seed = 122435;
        rng(seed);
        l0 = round(durpsp/tstep);

        if ~simpleTest
            dur = 1000; %ms
            xE = rand(round(rateE*dur*2),1);
            xI = rand(round(rateI*dur*2),1);
            tE = zeros(round(rateE*dur*2),1);
            tI = zeros(round(rateI*dur*2),1);
            pfE = randi(nE-ignorefE,round(rateE*dur*2),1);
            pfI = randi(nI-ignorefI,round(rateI*dur*2),1);
            t = 0;
            it = 0;
            if rateE > 0
                while t < dur
                    it = it + 1;
                    tE(it) = t;
                    t = t - log(xE(it))/rateE;
                end
                para.f_E = fE(pfE(1:it));
                nfE = it;
                para.tE = round(tE(1:it)'/tstep)*tstep;
            else
                para.tE = [];
                para.f_E = [];
            end
            t = 0;
            it = 0;
            if rateI > 0
                while t < dur
                    it = it +1;
                    tI(it) = t;
                    t = t - log(xI(it))/rateI;
                end
                para.f_I = fI(pfI(1:it));
                nfI = it;
                para.tI = round(tI(1:it)'/tstep)*tstep;
            else
                para.tI = [];
                para.f_I = [];
            end
            fname = ['MultipleInputs',num2str(dur),'-E',num2str(rateE*1e3),'-I',num2str(rateI*1e3),'-S',num2str(seed),'-',name,'-',num2str(i),'th-v'];
        else
            pfE = nE-ignorefE;
            pfI = nI-ignorefI;
            nfE = length(pfE);
            nfI = length(pfI);
            para.f_E = fE(pfE);
            para.f_I = fI(pfI);
            para.tE = 0;
            para.tI = 10;
            dur = dur0 + para.tI;
            fname = ['simpleTest',num2str(dur),'-E',num2str(para.tE),'-I',num2str(para.tI),'-',name,'-',num2str(i),'th-v'];
        end
        if test && exist([fname,'.mat'],'file')
            load([fname,'.mat'],'tmp');
        else
            tic;
            para.vtime = -1;
            switch model
                case 'HH'
                    tmp = silico(name,para.vRest(i),para,bool,tstep,dur,i,false);
                case 'EIF'
                    tmp = silico(name,para.vRest(i),para,bool,tstep,dur,i,false);
                case 'IF'
                    tmp = silico(name,para.vRest(i),para,bool,tstep,dur,i,false);
            end
            save([fname,'.mat'],'tmp');
            toc;
        end
        v = squeeze(tmp(1,:,:));
        vTargetE = v(round(para.tE/tstep)+1);
        vTargetI = v(round(para.tI/tstep)+1);
        t = 0:tstep:dur;
        nt = length(t);
        %ftE = floor(para.tE/tstep)+1;
        %ctE = ceil(para.tE/tstep)+1;
        %vTargetE = v(ftE) + (para.tE/tstep+1-ftE).*(v(ctE) - v(ftE));
        %ftI = floor(para.tI/tstep)+1;
        %ctI = ceil(para.tI/tstep)+1;
        %vTargetI = v(ftI) + (para.tI/tstep+1-ftI).*(v(ctI) - v(ftI));
        gE = zeros(size(v));
        gI = gE;

        vE = zeros(nt0,nfE);
        tEs = zeros(nfE,1);
        tEe = zeros(nfE,1);
        tEl = zeros(nfE,1);
        tEeadd = zeros(nfE,1);
        tEladd = zeros(nfE,1);
        vI = zeros(nt0,nfI);
        tIs = zeros(nfI,1);
        tIe = zeros(nfI,1);
        tIl = zeros(nfI,1);
        tIeadd = zeros(nfI,1);
        tIladd = zeros(nfI,1);
        vpred = zeros(nt,1);
        vadd = zeros(nt,1);
        for iE = 1:nfE
            tEs(iE) = round(para.tE(iE)/tstep)+1;
            tEe(iE) = min(tEs(iE)+l0,nt);
            tEl(iE) = tEe(iE)-tEs(iE)+1;
            tEeadd(iE) = min(tEs(iE)-1+nt0,nt);
            tEladd(iE) = tEeadd(iE)-tEs(iE)+1;
            tpick = tEs(iE):tEeadd(iE);
            gE(tpick) = gE(tpick) + G_E(0,t0(1:tEladd(iE)),fE(pfE(iE)),para.tau_er,para.tau_e,para.gConE);
        end
        for iI = 1:nfI
            tIs(iI) = round(para.tI(iI)/tstep)+1;
            tIe(iI) = min(tIs(iI)+l0,nt);
            tIl(iI) = tIe(iI)-tIs(iI)+1;
            tIeadd(iI) = min(tIs(iI)-1+nt0,nt);
            tIladd(iI) = tIeadd(iI)-tIs(iI)+1;
            tpick = tIs(iI):tIeadd(iI);
            gI(tpick) = gI(tpick) + G_E(0,t0(1:tIladd(iI)),fI(pfI(iI)),para.tau_ir,para.tau_i,para.gConI);
        end
        disp('linear add:');
        tic; 
        vadd = vadd + para.vRest(i);
        iI = 1;
        finishI = false;
        for iE = 1:nfE
            while para.tI(iI) < para.tE(iE)
                vI(:,iI) = interpPSP(sIPSP0(:,pfI(iI),:),vadd(tIs(iI)),vRange);
                tpick = tIs(iI):tIeadd(iI);
                vadd(tpick) = vadd(tpick) + vI(1:tIladd(iI),iI);
                if iI == nfI
                    finishI = true;
                    break;
                else
                    iI = iI + 1;
                end
            end
            vE(:,iE) = interpPSP(sEPSP0(:,pfE(iE),:),vadd(tEs(iE)),vRange);
            tpick = tEs(iE):tEeadd(iE);
            vadd(tpick) = vadd(tpick) + vE(1:tEladd(iE),iE);
        end
        if ~finishI
            iI0 = iI;
            for iI = iI0:nfI
                vI(:,iI) = interpPSP(sIPSP0(:,pfI(iI),:),vadd(tIs(iI)),vRange);
                tpick = tIs(iI):tIeadd(iI);
                vadd(tpick) = vadd(tpick) + vI(1:tIladd(iI),iI);
            end
        end
            %vpred = interpPSPE + interpPSPI;
        toc;
        disp('bilinear correction:');
        tic;
        bEI = cell(nfE,nfI);
        tEI = zeros(2,nfE,nfI);  
        bEE = cell(nfE,nfE);
        tEE = zeros(2,nfE,nfI);
        bIE = cell(nfI,nfE);
        tIE = zeros(2,nfI,nfE);
        bII = cell(nfI,nfI);
        tII = zeros(2,nfI,nfI);

        vpred = vpred + para.vRest(i);

        finishI = false;
        iI = 1;
%         debug = true;
        debug = false;
        if debug
            dbstop if error
            dbstop at 1259
            hDebug = figure;
            xDebug = [0,300];
            subplot(2,1,1);
            hold on
            plot(t,v,'k');
            plot(t,vadd,':m');
        end
        rule = @(iiE,iiI,iE,iI) (iiI == 1 && iE == 1) || (iiE == 1 && iI ==1);
        %rule = @(iiE,iiI,iE,iI) iI == 10;
        for iE = 1:nfE
            while para.tI(iI) < para.tE(iE) && ~finishI
                vI(:,iI) = interpPSP(sIPSP0(:,pfI(iI),:),vpred(tIs(iI)),vRange);
                tpick = tIs(iI):tIeadd(iI);
                vpred(tpick) = vpred(tpick) + vI(1:tIladd(iI),iI);
                if testEI
                    for iiE = 1: iE-1
                        if round((para.tI(iI) - para.tE(iiE))/tstep) < l0
                            if debug && rule(iiE,0,0,iI)
                                hDebug = debugPreplot(hDebug,t,vpred);
                            end
                            [vpred, bEI{iiE,iI},tEI(:,iiE,iI)] = pred(vpred,iiE,iI,squeeze(sEPSP(:,pfE(iiE),:,:)),vI,...
                                                                      tEs,tEe,tEl,tIs,kVEI,...
                                                                      vRange,nv0,ndt,idtRange);
                            if debug && rule(iiE,0,0,iI)
                                hDebug = plotDebug(hDebug,t,v,vpred,iiE,iI,...
                                                   bEI,tEI,pfE,pfI,...
                                                   tEs,tEe,tIs,tIe,tIl,'EI',xDebug,...
                                                   vE,vI);
                                rule = @(iiE,iiI,iE,iI) true;
                            end
                        end
                    end
                end
                if testII
                    for iiI = 1: iI-1
                        if round((para.tI(iI) - para.tI(iiI))/tstep) < l0
                            if debug && rule(0,iiI,0,iI)
                                hDebug = debugPreplot(hDebug,t,vpred);
                            end
                            [vpred, bII{iiI,iI},tII(:,iiI,iI)] = pred(vpred,iiI,iI,squeeze(sIPSP(:,pfI(iiI),:,:)),vI,...
                                                                      tIs,tIe,tIl,tIs,kVII,...
                                                                      vRange,nv0,ndt,idtRange);
                            if debug && rule(0,iiI,0,iI)
                                hDebug = plotDebug(hDebug,t,v,vpred,iiI,iI,...
                                                   bII,tII,pfI,pfI,...
                                                   tIs,tIe,tIs,tIe,tIl,'II',xDebug,...
                                                   vI,vI);
                                rule = @(iiE,iiI,iE,iI) true;
                            end
                        end
                    end
                end
                if iI == nfI
                    finishI = true;
                    break;
                else
                    iI = iI + 1;
                end
            end

            vE(:,iE) = interpPSP(sEPSP0(:,pfE(iE),:),vpred(tEs(iE)),vRange);
            tpick = tEs(iE):tEeadd(iE);
            vpred(tpick) = vpred(tpick) + vE(1:tEladd(iE),iE);
            if testIE
                for iiI = 1:iI-1
                    if round((para.tE(iE) - para.tI(iiI))/tstep) < l0 && round((para.tE(iE) - para.tI(iiI))/tstep) > 0
                        if debug && rule(0,iiI,iE,0)
                            hDebug = debugPreplot(hDebug,t,vpred);
                        end
                        [vpred, bIE{iiI,iE},tIE(:,iiI,iE)] = pred(vpred,iiI,iE,squeeze(sIPSP(:,pfI(iiI),:,:)),vE,...
                                                                  tIs,tIe,tIl,tEs,kVIE,...
                                                                  vRange,nv0,ndt,idtRange);
                        if debug && rule(0,iiI,iE,0)
                            hDebug = plotDebug(hDebug,t,v,vpred,iiI,iE,...
                                               bIE,tIE,pfI,pfE,...
                                               tIs,tIe,tEs,tEe,tEl,'IE',xDebug,...
                                               vI,vE);
                            rule = @(iiE,iiI,iE,iI) true;
                        end
                    end
                end
            end
            if testEE
                for iiE = 1:iE-1
                    if round((para.tE(iE) - para.tE(iiE))/tstep) < l0
                        if debug && rule(iiE,0,iE,0)
                            hDebug = debugPreplot(hDebug,t,vpred);
                        end
                        [vpred, bEE{iiE,iE},tEE(:,iiE,iE)] = pred(vpred,iiE,iE,squeeze(sEPSP(:,pfE(iiE),:,:)),vE,...
                                                                  tEs,tEe,tEl,tEs,kVEE,...
                                                                  vRange,nv0,ndt,idtRange);
                        if debug && rule(iiE,0,iE,0)
                            hDebug = plotDebug(hDebug,t,v,vpred,iiE,iE,...
                                               bEE,tEE,pfE,pfE,...
                                               tEs,tEe,tEs,tEe,tEl,'EE',xDebug,...
                                               vE,vE);
                            rule = @(iiE,iiI,iE,iI) true;
                        end
                    end
                end
            end
        end
        if ~finishI
            iI0 = iI;
            for iI = iI0:nfI

                vI(:,iI) = interpPSP(sIPSP0(:,pfI(iI),:),vpred(tIs(iI)),vRange);
                tpick = tIs(iI):tIeadd(iI);
                vpred(tpick) = vpred(tpick) + vI(1:tIladd(iI),iI);
                if testEI
                    for iiE = 1: nfE
                        if round((para.tI(iI) - para.tE(iiE))/tstep) < l0 
                            if debug && rule(iiE,0,0,iI)
                                hDebug = debugPreplot(hDebug,t,vpred);
                            end
                            [vpred,bEI{iiE,iI},tEI(:,iiE,iI)] = pred(vpred,iiE,iI,squeeze(sEPSP(:,pfE(iiE),:,:)),vI,...
                                                                     tEs,tEe,tEl,tIs,kVEI,...
                                                                     vRange,nv0,ndt,idtRange);
                            if debug && rule(iiE,0,0,iI)
                                hDebug = plotDebug(hDebug,t,v,vpred,iiE,iI,...
                                                   bEI,tEI,pfE,pfI,...
                                                   tEs,tEe,tIs,tIe,tIl,'EI',xDebug,...
                                                   vE,vI);
                                rule = @(iiE,iiI,iE,iI) true;
                            end
                        end
                    end
                end
                if testII
                    for iiI = 1: iI-1
                        if round((para.tI(iI) - para.tI(iiI))/tstep) < l0
                            if debug && rule(0,iiI,0,iI)
                                hDebug = debugPreplot(hDebug,t,vpred);
                            end
                            [vpred,bII{iiI,iI},tII(:,iiI,iI)] = pred(vpred,iiI,iI,squeeze(sIPSP(:,pfI(iiI),:,:)),vI,...
                                                                     tIs,tIe,tIl,tIs,kVII,...
                                                                     vRange,nv0,ndt,idtRange);
                            if debug && rule(0,iiI,0,iI)
                                hDebug = plotDebug(hDebug,t,v,vpred,iiI,iI,...
                                                   bII,tII,pfI,pfI,...
                                                   tIs,tIe,tIs,tIe,tIl,'II',xDebug,...
                                                   vI,vI);
                                rule = @(iiE,iiI,iE,iI) true;
                            end
                        end
                    end
                end
            end
        end
        toc;
        hM = figure;
        ax1 = subplot(2,1,1);
        hold on
        targetVersion = '2016b';
        currentVersion = version('-release');
        compatible = currentVersion>targetVersion;
        vmean = mean(v);
        if sum(compatible) > 0 || strcmp(currentVersion,targetVersion)
            yyaxis right
            hold on            
            plot(t,gE*(para.vE-vmean),':r');
            plot(t,gI*(vmean-para.vI),':b');
            plot(t,gE.*(para.vE-v),'--r');
            plot(t,gI.*(v-para.vI),'--b');
            if strcmp(model,'HH')
                iNa = squeeze(tmp(2,:,:)).^3.*squeeze(tmp(4,:,:)).*para.gNa(i).*(v-para.vNa);
                iKd = squeeze(tmp(3,:,:)).^4.*para.gK(i).*(v-para.vK);
                plot(t,abs(iNa),'-r');
                plot(t,abs(iKd),'-b');
            end
            yyaxis left
            hold on
            plot(t,v);
        else
            [ax,h1,h2] = plotyy(t,v,t,gE*(para.vE-vmean));
            h2.LineStyle = ':';
            h2.Color = 'r';
            axes(ax(2));
            hold(ax(2),'on');
            plot(t,gI*(vmean-para.vI),':b');
            plot(t,gE.*(para.vE-v),'--r');
            plot(t,gI.*(v-para.vI),'--b');
            if strcmp(model,'HH')
                iNa = squeeze(tmp(2,:,:)).^3.*squeeze(tmp(4,:,:)).*para.gNa(i).*(v-para.vNa);
                iKd = squeeze(tmp(3,:,:)).^4.*para.gK(i).*(v-para.vK);
                plot(t,abs(iNa),'-r');
                plot(t,abs(iKd),'-b');
            end
            ylim(ax(2),[0,inf]);
            axes(ax(1));
            hold(ax(1),'on');
            ax(1).Color = 'none';
        end
        plot(t,vadd,'--m');
        plot(t,vpred,':m');
        sl = (max(v)-min(v));
        y = [min(v)-sl*1.2,max(v) + sl*1.2];
        FS = 10;
        plot(ones(2,1)*para.tE,[y(1)*ones(1,nfE);vTargetE],':r');
        plot(para.tE,y(1)*ones(1,nfE),'sr');
        plot(ones(2,1)*para.tI,[y(1)*ones(1,nfI);vTargetI],':b');
        plot(para.tI,y(1)*ones(1,nfI),'ob');
        xl = [t(1),t(end)];
        if poiend
            selectE = tEs + l0 < nt;
            sfE = sum(selectE);
            selectI = tIs + l0 < nt;
            sfI = sum(selectI);
            plot(ones(2,1)*para.tE(selectE)+durpsp,[v(tEs(selectE)+l0)';y(2)*ones(1,sfE)],':r');
            plot(ones(2,1)*para.tI(selectI)+durpsp,[v(tIs(selectI)+l0)';y(2)*ones(1,sfI)],':b');
        end
        if poinumber
            for iE = 1:nfE
                if para.tE(iE) > xl(1) && para.tE(iE) < xl(2)
                    text(para.tE(iE),vTargetE(iE)-sl*0.1,num2str(iE),'Color','r','FontSize',FS);
                    text(para.tE(iE),vTargetE(iE)-sl*0.2,num2str(pfE(iE)),'Color','r','FontSize',FS);
                end
                if poiend
                    if para.tE(iE)+durpsp > xl(1) && para.tE(iE)+durpsp < xl(2)
                        text(para.tE(iE)+durpsp,v(tEs(iE)+l0)+sl*0.1,num2str(iE),'Color','r','FontSize',FS);
                        text(para.tE(iE)+durpsp,v(tEs(iE)+l0)+sl*0.2,num2str(pfE(iE)),'Color','r','FontSize',FS);
                    end
                end
            end
            for iI = 1:nfI
                if para.tI(iI) > xl(1) && para.tI(iI) < xl(2)
                    text(para.tI(iI),vTargetI(iI)-sl*0.1,num2str(iI),'Color','b','FontSize',FS);
                    text(para.tI(iI),vTargetI(iI)-sl*0.2,num2str(pfI(iI)),'Color','b','FontSize',FS);
                end
                if poiend
                    if para.tI(iI)+durpsp > xl(1) && para.tI(iI)+durpsp < xl(2)
                        text(para.tI(iI)+durpsp,v(tIs(iI)+l0)+sl*0.1,num2str(iI),'Color','b','FontSize',FS);
                        text(para.tI(iI)+durpsp,v(tIs(iI)+l0)+sl*0.2,num2str(pfI(iI)),'Color','b','FontSize',FS);
                    end
                end
            end
        end
        plot(t(1)*ones(1,nv0),vRange,'dk');
        ylim(y);
        xlim(xl);
        xlabel('ms');
        ylabel('mV');
        ax3 = subplot(2,2,3);
        hold on
        v = v';
        plot([0,3],zeros(2,1),':k');
        errorbar(1,sign(sum(vadd -v>0)*2-1) *mean(abs(vadd-v)),  std(abs(vadd-v)),'*m');
        e=errorbar(2,sign(sum(vpred-v>0)*2-1)*mean(abs(vpred-v)),std(abs(vpred-v)),'*b');
        set(ax3,'XTick',[1,2]);
        set(ax3,'XTickLabel',{'linear';'bilinear'});
        legend([e],{'kV'});
        ylabel('error mV');
        xlim([0.5,2.5]);
        ax4 = subplot(2,2,4);
%         xl = [t(min(tEs(round(nfE/2)),tIs(round(nfI/2)))),t(nt)];
        ix = min(23,nfE);
        ixl = [tEs(ix),tEe(ix)];
%         ixl = [min(tEs(round(nfE/2)),tIs(round(nfI/2))),min(tEe(round(nfE*3/4)),tIe(round(nfI*3/4)))];
        xl = t(ixl);
        hs = findobj(ax1,'-not','Type','Axes','-and','-not','Type','Text');
        copyobj(hs,ax4);
        hold on
        sl = max(v(ixl(1):ixl(2)))-min(v(ixl(1):ixl(2)));
        y2 = [min(v(ixl(1):ixl(2)))-sl*1.2,max(v(ixl(1):ixl(2))) + sl*1.2];
        for iE = 1:nfE
            if para.tE(iE) > xl(1) && para.tE(iE) < xl(2)
                text(para.tE(iE),vTargetE(iE)-sl*0.1,num2str(iE),'Color','r','FontSize',FS);
                text(para.tE(iE),vTargetE(iE)-sl*0.2,num2str(pfE(iE)),'Color','r','FontSize',FS);
            end
            if para.tE(iE)+durpsp > xl(1) && para.tE(iE)+durpsp < xl(2)
                text(para.tE(iE)+durpsp,v(tEs(iE)+l0)+sl*0.1,num2str(iE),'Color','r','FontSize',FS);
                text(para.tE(iE)+durpsp,v(tEs(iE)+l0)+sl*0.2,num2str(pfE(iE)),'Color','r','FontSize',FS);
            end
        end
        for iI = 1:nfI
            if para.tI(iI) > xl(1) && para.tI(iI) < xl(2)
                text(para.tI(iI),vTargetI(iI)-sl*0.1,num2str(iI),'Color','b','FontSize',FS);
                text(para.tI(iI),vTargetI(iI)-sl*0.2,num2str(pfI(iI)),'Color','b','FontSize',FS);
            end
            if para.tI(iI)+durpsp > xl(1) && para.tI(iI)+durpsp < xl(2)
                text(para.tI(iI)+durpsp,v(tIs(iI)+l0)+sl*0.1,num2str(iI),'Color','b','FontSize',FS);
                text(para.tI(iI)+durpsp,v(tIs(iI)+l0)+sl*0.2,num2str(pfI(iI)),'Color','b','FontSize',FS);
            end
        end
        if ~poiend
            selectE = (para.tE + durpsp <= xl(2) ) & (para.tE + durpsp >= xl(1));
            sfE = sum(selectE);
            selectI = (para.tI + durpsp <= xl(2) ) & (para.tI + durpsp >= xl(1));
            sfI = sum(selectI);
            if sfE>0
                plot(ones(2,1)*para.tE(selectE)+durpsp,[v(tEs(selectE)+l0)';y(2)*ones(1,sfE)],':r');
            end
            if sfI>0
                plot(ones(2,1)*para.tI(selectI)+durpsp,[v(tIs(selectI)+l0)';y(2)*ones(1,sfI)],':b');
            end
        end
        xlabel('ms');
        ylabel('mV');
        xlim(xl);
        
        ylim([y2(1),y2(2)]);
        printpic(hM,diri,fname,picformat,printDriver,dpi,pos0);
    end
end
function [EPSP,IPSP,E_tmax,I_tmax,vleakage] = sPSP_check(silico,tstep,vRange,fIRange,fERange,para,bool,name,dur,i,v0id,dtRange)
    ndt = length(dtRange);
    nv0 = length(vRange);
    nE = length(fERange);
    nI = length(fIRange);
    nt = round(dur/tstep)+1;
    EPSP = zeros(nt,nE,ndt,nv0);
    E_tmax =  zeros(nE,ndt,nv0);
    IPSP = zeros(nt,nI,ndt,nv0);
    I_tmax = zeros(nI,ndt,nv0);
    vleakage = zeros(nt,ndt,nv0);

    v0 = vRange(v0id);
    para.f_E = 0;
    para.tE = para.f_E;
    para.f_I = 0;
    para.tI = para.f_I;
    for iv0 = 1:nv0
        para.newv = vRange(iv0);
        parfor idt = 1:ndt
            param = para;
            param.vtime = dtRange(idt);
            leakage = silico(name,v0,param,bool,tstep,dur,i,false);
            vleakage(:,idt,iv0) = squeeze(leakage(1,:,:));
        end
    end
    disp('leakage complete');
    for iv0 = 1:nv0
        para.newv = vRange(iv0);

        v0 = vRange(v0id)*ones(nE,1);
        para.f_E = fERange;
        para.tE = zeros(nE,1);
        para.f_I = para.tE;
        para.tI = para.tE;

        parfor idt = 1:ndt
            param = para;
            param.vtime = dtRange(idt);
            tmpE = silico(name,v0,param,bool,tstep,dur,i,false);
            tmpEv = squeeze(tmpE(1,:,:)) - repmat(vleakage(:,idt,iv0),[1,nE]);
            [~,E_tmax(:,idt,iv0)] = max(tmpEv);
            EPSP(:,:,idt,iv0) = tmpEv;
        end
        for idt = 1:ndt
            for it = 1:nt
                for iF = 1:nE
                    assert(EPSP(it,iF,idt,iv0)>=0);
                end
            end
        end

        v0 = vRange(v0id)*ones(nI,1);
        para.f_I = fIRange;
        para.tI = zeros(nI,1);
        para.f_E = para.tI;
        para.tE = para.tI;
         
        parfor idt = 1:ndt
            param = para;
            param.vtime = dtRange(idt);
            tmpI = silico(name,v0,param,bool,tstep,dur,i,false);
            tmpIv = squeeze(tmpI(1,:,:)) - repmat(vleakage(:,idt,iv0),[1,nI]);
            [~,I_tmax(:,idt,iv0)] = min(tmpIv);
            IPSP(:,:,idt,iv0) = tmpIv;
        end
        for idt = 1:ndt
            for it = 1:nt
                for iF = 1:nI
                    assert(IPSP(it,iF,idt,iv0)<=0);
                end
            end
        end
    end
    disp('sPSP complete');
end
function [h,EPSP,IPSP,vleakage] = sPSP0_check(silico,tstep,vRange,fIRange,fERange,para,bool,name,dur,i,v0id)
    nv0 = length(vRange);
    nE = length(fERange);
    nI = length(fIRange);
    nt = round(dur/tstep)+1;
    EPSP = zeros(nt,nE,nv0);
    IPSP = zeros(nt,nI,nv0);
    vleakage = zeros(nt,nv0);

    para.f_E = 0;
    para.tE = para.f_E;
    para.f_I = 0;
    para.tI = para.f_I;
    para.vtime = -1;
    parfor iv0 = 1:nv0
        param = para;
        v0 = vRange(iv0);
        leakage = silico(name,v0,param,bool,tstep,dur,i,false);
        vleakage(:,iv0) = squeeze(leakage(1,:,:));
    end
    disp('leakage0 complete');
    parfor iv0 = 1:nv0
       param = para;

       v0 = vRange(iv0)*ones(nE,1);
       param.f_E = fERange;
       param.tE = zeros(nE,1);
       param.f_I = param.tE;
       param.tI = param.tE;

       tmpE = silico(name,v0,param,bool,tstep,dur,i,false);
       EPSP(:,:,iv0) = squeeze(tmpE(1,:,:)) - repmat(vleakage(:,iv0),[1,nE]);

       v0 = vRange(iv0)*ones(nI,1);
       param.f_I = fIRange;
       param.tI = zeros(nI,1);
       param.f_E = param.tI;
       param.tE = param.tI;
        
       tmpI = silico(name,v0,param,bool,tstep,dur,i,false);
       IPSP(:,:,iv0) = squeeze(tmpI(1,:,:)) - repmat(vleakage(:,iv0),[1,nI]);
    end
    disp('sPSP0 complete');
    h = figure;
    s = [0.1,1.0];
    v = 1.0;
    t = linspace(0,dur,nt);
    subplot(2,2,2);
    plotSinglePSP(0,s,v,fERange,t,EPSP,nE,nv0,v0id);
    xlabel('f');
    ylabel('t');
    zlabel('EPSP (mV)');
    subplot(2,2,4);
    plotSinglePSP(2/3,s,v,fIRange,t,IPSP,nI,nv0,v0id);
    xlabel('f');
    ylabel('t');
    zlabel('IPSP (mV)');
    c = zeros(nv0,3);
    c(:,3) = linspace(0.7,0,nv0)';
    c = hsv2rgb(c);
    t = 0:tstep:dur;
    subplot(1,2,1);
    hold on
    for iv0 = 1:nv0
        plot(t,squeeze(vleakage),'Color',c(iv0,:));
    end
    xlabel('t');
    ylabel('Membrane Potential (mV)');
    title('leakage');
end
function h = plotsPSP0(EPSP,IPSP,vleakage,fERange,fIRange,nv0,dur,nt,v0id)
    nE = length(fERange);
    nI = length(fIRange);
    h = figure;
    s = [0.1,1.0];
    v = 1.0;
    t = linspace(0,dur,nt);
    subplot(2,2,2);
    plotSinglePSP(0,s,v,fERange,t,EPSP,nE,nv0,v0id);
    xlabel('f');
    ylabel('t');
    zlabel('EPSP (mV)');
    subplot(2,2,4);
    plotSinglePSP(2/3,s,v,fIRange,t,IPSP,nI,nv0,v0id);
    xlabel('f');
    ylabel('t');
    zlabel('IPSP (mV)');
    c = zeros(nv0,3);
    c(:,3) = linspace(0.7,0,nv0)';
    c = hsv2rgb(c);
    subplot(1,2,1);
    hold on
    for iv0 = 1:nv0
        plot(t,squeeze(vleakage),'Color',c(iv0,:));
    end
    xlabel('t');
    ylabel('Membrane Potential (mV)');
    title('leakage');
end
function plotSinglePSP(h,s,v,x,y,z,m,n,j0)
    hold on
    grid on
    nstep = length(y);
    c = zeros(n,3);
    c(:,1) = h;
    c(:,2) = linspace(s(1),s(2),n)';
    c(:,3) = v;
    c = hsv2rgb(c);
    for i = 1:m
        xx = zeros(nstep,1)+x(i);
        for j = 1:n
            if j == j0
                plot3(xx,y,z(:,i,j),'k');
            else
                plot3(xx,y,z(:,i,j),'Color',c(j,:));
            end
        end
    end
    [mz, ind] = max(abs(z(:)));
    campos([0.1*max(x),0.2*max(y),0.8*sign(z(ind))*mz]);
end
function [v,v1,v2,ind] = interpPSP(sPSP,vTarget,vRange)
    nv0 = length(vRange); 
    tmp = (vTarget-vRange)<=0;
    ind =find(tmp,1,'first');
    if isempty(ind)
        v1 = sPSP(:,:,nv0-1);
        v2 = sPSP(:,:,nv0);
        v = v1 + (vTarget-vRange(nv0-1))/(vRange(nv0)-vRange(nv0-1))*(v2-v1);
    else
        if (ind == 1)
            v1 = sPSP(:,:,1);
            v2 = sPSP(:,:,2);
            v = v2 + (vTarget-vRange(2))/(vRange(1)-vRange(2))*(v1-v2);
        else
            v1 = sPSP(:,:,ind-1);
            v2 = sPSP(:,:,ind);
            v = v1 + (vTarget-vRange(ind-1))/(vRange(ind)-vRange(ind-1))*(v2-v1);
        end
    end
end
function tmp = interpPSPdt(vec,i,tar,ref,l,eov)
    nref = length(ref);
    if tar >= ref(nref)
        i = nref-1;
        j = nref;
    else
        while ref(i+1) <  tar
            i = i + 1;  
        end
        j = i+1;
    end
    r = (tar-ref(i))/(ref(j)-ref(i));
    base = vec(tar+(0:l-1),i);
    l0 = eov-ref(j)+1;
    if l0 < l
        tmp = base + r*([vec(ref(j):eov,j);zeros(l-l0,1)]-base);
    else
        tmp = base + r*(vec(ref(j)+(0:l-1),j)-base);
    end
end
function [kV,pV,vaddV,vDoubletV,EI,h0] = doubleCheck(silico,name,para,v0,bool,tstep,iidt,dur,i,vleakage,sPSP1,sPSP2,vRange,pp,idt,ndt,dtRange)
    global nvplot iv0Case
    n1 = size(sPSP1,2);
    n2 = size(sPSP2,2);
    nt = round(dur/tstep)+1;
    nv0 = length(vRange);
    i2 = n2;
    i1 = n1;

    kV = zeros(nt,ndt,nv0);
    vDoubletV = zeros(nt,n1*n2,nv0);
    vaddV = zeros(n1*n2,nt,nv0);
    v1V = zeros(n1*n2,nt,nv0);
    v2V = zeros(n1*n2,nt,nv0);
    pV = zeros(n1*n2,nt,nv0);
    para.vtime = (iidt-1)*tstep;
    dtOI = iidt:nt;
    dtl = length(dtOI);
    parfor iv0 = 1:nv0
       param = para;
       param.newv = vRange(iv0);
       tmp = silico(name,v0,param,bool,tstep,dur,i,false);
       vDoubletV(:,:,iv0) = squeeze(tmp(1,:,:));
       vDoubletV(:,:,iv0) = vDoubletV(:,:,iv0) - repmat(vleakage(:,idt,iv0),[1,n1*n2]);
    end
    for iv0 = 1:nv0
       for iF = 1:n1
           range = (iF-1)*n2+(1:n2);
           v1V(range,dtOI,iv0) = (sPSP1(dtOI,iF,idt,iv0)*ones(1,n2))';
           v2V(range,dtOI,iv0) = [sPSP2(1:dtl,:,iv0)]';
       end
       vaddV(:,dtOI,iv0) = v1V(:,dtOI,iv0) + v2V(:,dtOI,iv0);
       parfor it = iidt:nt
       %for it = 1:nt
           [kV(it,1,iv0),pV(:,it,iv0)] = p_fit110k0(v1V(:,it,iv0),v2V(:,it,iv0),vDoubletV(it,:,iv0)');
       end
    end

    EI = v1V.*v2V;
    h0 = 0;
    if pp 
        h0 = figure;
        t = 0:tstep:dur;
        range = (i1-1)*n2+(1:n2);
        for iv =1:nvplot
            iv0 = iv0Case(iv); 
            subplot(2,nvplot,iv);
            title([num2str(iv),'th v']);
            hold on
            pick = iidt:nt;
            plot(t(pick),v1V(range(1),pick,iv0),':r');
            plot(t(pick),v2V(range(i2),pick,iv0),':b');
            plot(t(pick),vaddV(range(i2),pick,iv0),':k');
            plot(t(pick),pV(range(i2),pick,iv0),'--m');
            [ax,h1,h2] = plotyy(t(pick),vDoubletV(pick,range(i2),iv0),t(pick),kV(pick,1,iv0));
            h1.Color = 'k';
            h2.LineStyle = ':';
            h2.Color = 'g';
            subplot(2,nvplot,nvplot+iv);
            [ax,h1,h2] = plotyy(t(pick),h(pick,range(i2),iv0),t(pick),[m(pick,range(i2),iv0),n(pick,range(i2),iv0)]);
            h1.Color = 'g';
            h2(1).Color = 'r';
            h2(2).Color = 'b';
        end
    end
    idtRange = round(dtRange/tstep)+1;
    vtmp = zeros(nt,n1*n2,nv0);
    for jdt = idt+1:ndt
        jjdt = idtRange(jdt);
        kkdt = (iidt+jjdt-1);
        dtOI = kkdt:nt;
        dtl = length(dtOI);
        para.vtime = kkdt*tstep;
        parfor iv0 = 1:nv0
            param = para;
            param.newv = vRange(iv0);
            tmp = silico(name,v0,param,bool,tstep,dur,i,false);
            vtmp(:,:,iv0) = squeeze(tmp(1,:,:));
        end
        for iv0 = 1:nv0
            tmpv = vleakage(1:dtl,1,iv0);
            vtmp(dtOI,:,iv0) = vtmp(dtOI,:,iv0) - repmat(tmpv,[1,n1*n2]);
            for iF = 1:n1
                range = (iF-1)*n2+(1:n2);
                tmp = interpPSPdt(squeeze(sPSP1(:,iF,:,iv0)),idt,kkdt,idtRange,dtl,nt);
                v1V(range,dtOI,iv0) = (tmp*ones(1,n2))';
                v2V(range,dtOI,iv0) = [sPSP2(1:dtl,:,iv0)]';
            end
            parfor it = kkdt:nt
            %for it = 1:nt
                [kV(it,jdt,iv0),~] = p_fit110k0(v1V(:,it,iv0),v2V(:,it,iv0),vtmp(it,:,iv0)');
            end
        end
    end
end
function h = drawExample(h,pv0,vadd0,vDoublet0,vleakage0,t,iidt,tp0,k0,n,vRange,vRest,v120,xl,titl)
    global nvplot iv0Case;
    nv0 = length(vRange);
    inv0 = 0;
    for iv0 = iv0Case;
        vleakage = repmat(vleakage0(:,iv0),[1,n]);
        k = k0(:,iv0);
        v12 = v120(:,:,iv0);
        vDoublet = vDoublet0(:,:,iv0) + vleakage - vRest;
        pv = pv0(:,:,iv0) + vleakage' - vRest;
        vadd = vadd0(:,:,iv0) + vleakage' - vRest;
        nt = length(t);
        figure(h);
        tp = min(iidt + tp0,nt);
        subplot(2,2,1);
        hold on
        for iF=1:n
            p = plot(t,vDoublet(:,iF));
            plot(t,pv(iF,:),':','Color',p.Color,'LineWidth',1.5);
            plot(t,vadd(iF,:),'--','Color',p.Color);
        end
        plot(t(tp),vDoublet(tp,:),'*');
        plot(ones(2,1)*t(iidt),ylim,':k');
        xlabel('time (ms)');
        ylabel('PSP (mv)');
        title([titl]);

        subplot(2,2,3);
        hold on
        plot(vDoublet(tp,:)',pv(:,tp),'o');
        plot(vDoublet(tp,:)',vadd(:,tp),'s');
        x = xlim;
        y = x;
        plot(x,y,':k');
        xlabel('PSP');
        ylabel('predict');

        subplot(nvplot,2,2*inv0+2)
        hold on
        plot(v12(:,tp),(vDoublet(tp,:)'-vadd(:,tp)),'o');
        xx = sort(v12(:,tp));
        plot(xx,xx*k(tp),':k','LineWidth',1.5);
        ylabel('V_{SC}');
        xlabel(xl);
        title({['k = ',num2str(k(tp),'%.3g')]});
        inv0 = inv0 + 1;
    end
end
function printpic(h,dir,fname,picformat,printDriver,dpi,pos)
    if ~isempty(picformat)
        set(h,'PaperUnits', 'inches','PaperSize',pos(1:2));
        set(h,'PaperPosition',pos(3:6));
    %         set(h,'PaperType','usletter');
    %         set(h,'PaperOrientation','portrait');
        set(h,'Units','inches');
    %     set(h,'OuterPosition',pos(7:10));
    %     set(h,'Position',pos(3:6)+8);
        set(h,'Renderer','Painter');
    % %     set(h,'PaperPositionMode','auto');
        if ~strcmp(picformat,'fig')
            print(h,[dir,'/',fname,'.',picformat],printDriver,'-loose',dpi);
        end
        saveas(h,[dir,'/',fname,'.fig']);
%         close(h);
    end
end
function [vpred, b, t] = pred(vpred,i1,i2,sPSP1,v2,...
                              t1s,t1e,t1l,t2s,k,...
                              vRange,nv0,ndt,idtRange)
    lt = t1e(i1)-t2s(i2);
    assert(lt>0);
    tpick = t2s(i2)+(0:lt);
    b = zeros(length(tpick),2);
    [kVtmp,v1] = interpKV(k,sPSP1,vpred(t2s(i2)),t2s(i2)-t1s(i1)+1,vRange,nv0,idtRange,ndt,lt);
    b(:,1) = kVtmp.*v1.*v2(1:lt+1,i2);
    b(:,2) = kVtmp;
    vpred(tpick) = vpred(tpick) + b(:,1);
    t(1) = tpick(1); t(2) = tpick(end);
end

function [kVtmp, vtmp] = interpKV(k,sPSP,vTarget,tTarget,v,nv0,idtRange,ndt,lt)
    fv = find(vTarget-v>=0, 1,'last');
    fdt = find(tTarget-idtRange>=0, 1,'last');
    assert(~isempty(fdt));
    assert(fdt < ndt);

    t0 = (idtRange(fdt)) + (0:lt);
    t1 = (idtRange(fdt+1)) + (0:lt);
    rt = (tTarget-idtRange(fdt))/(idtRange(fdt+1)-idtRange(fdt));
    if isempty(fv) % out of lower bound exterpolation
        k0 = k(t0,fdt,1);
        kVtmp = k0 + (vTarget-v(1))/(v(1)-v(2)) * (k0-k(t0,fdt,2)) + rt * (k(t1,fdt+1,1)-k0);
        sPSP0 = sPSP(t0,fdt,1);
        vtmp = sPSP0 + (vTarget-v(1))/(v(1)-v(2)) * (sPSP0-sPSP(t0,fdt,2)) + rt * (sPSP(t1,fdt+1,1)-sPSP0);
    else
        k0 = k(t0,fdt,fv);
        sPSP0 = sPSP(t0,fdt,fv);
        if fv == nv0  % out of upper bound exterpolation
            kVtmp = k0 + (vTarget-v(fv))/(v(fv)-v(fv-1)) * (k0-k(t0,fdt,fv-1));
            vtmp = sPSP0 + (vTarget-v(fv))/(v(fv)-v(fv-1)) * (sPSP0-sPSP(t0,fdt,fv-1));
        else
            kVtmp = k0 + (vTarget-v(fv))/(v(fv+1)-v(fv)) * (k(t0,fdt,fv+1)-k0);
            vtmp = sPSP0 + (vTarget-v(fv))/(v(fv+1)-v(fv)) * (sPSP(t0,fdt,fv+1)-sPSP0);
        end
        kVtmp = kVtmp + rt * (k(t1,fdt+1,fv)-k0);
        vtmp = vtmp + rt * (sPSP(t1,fdt+1,fv)-sPSP0);
    end
end
function h = debugPreplot(h,t,vpred)
    figure(h)
    subplot(2,1,1)
    hold on
    plot(t,vpred,'-b');
end
function h = plotDebug(h,t,v,vpred,i1,i2,...
                       b,tl,f1,f2,...
                       ts1,te1,ts2,te2,tl2,label,xl,...
                       v1,v2)
    switch label
        case 'EE'
            s = {'sr','or'};
            c = {'r','r'};
        case 'EI'
            s = {'sr','ob'};
            c = {'r','b'};
        case 'II'
            s = {'sb','ob'};
            c = {'b','b'};
        case 'IE'
            s = {'sb','or'};
            c = {'b','r'};
    end
    figure(h);
    ax1 = subplot(2,1,1);
    hold on
    plot(t,vpred,'--b');
    y=ylim;
    yl = y(2)-y(1);

    plot(ones(1,2)*t(ts1(i1)),[y(1),v(ts1(i1))],':','Color',c{1});
    plot(ones(1,2)*t(te1(i1)),[y(2),v(te1(i1))],':','Color',c{1});
    plot(t(ts1(i1)),v(ts1(i1)),s{1});
    plot(t(te1(i1)),v(te1(i1)),s{1});
    text(t(ts1(i1)),v(ts1(i1))-yl*0.1,num2str(i1),'FontSize',7,'Color',c{1});
    text(t(ts1(i1)),v(ts1(i1))-yl*0.2,num2str(f1(i1)),'FontSize',7,'Color',c{1});
    text(t(te1(i1)),v(te1(i1))+yl*0.1,num2str(i1),'FontSize',7,'Color',c{1});
    text(t(te1(i1)),v(te1(i1))+yl*0.2,num2str(f1(i1)),'FontSize',7,'Color',c{1});

    plot(ones(1,2)*t(ts2(i2)),[y(1),v(ts2(i2))],':','Color',c{2});
    plot(ones(1,2)*t(te2(i2)),[y(2),v(te2(i2))],':','Color',c{2});
    plot(t(ts2(i2)),v(ts2(i2)),s{2});
    plot(t(te2(i2)),v(te2(i2)),s{2});
    text(t(ts2(i2)),v(ts2(i2))-yl*0.1,num2str(i2),'FontSize',7,'Color',c{2});
    text(t(ts2(i2)),v(ts2(i2))-yl*0.2,num2str(f2(i2)),'FontSize',7,'Color',c{2});
    text(t(te2(i2)),v(te2(i2))+yl*0.1,num2str(i2),'FontSize',7,'Color',c{2});
    text(t(te2(i2)),v(te2(i2))+yl*0.2,num2str(f2(i2)),'FontSize',7,'Color',c{2});

    xlim(xl)
    title(label);

    ax2 = subplot(2,2,3);
    vl = te2(i2)-ts1(i1)+1;
    va = zeros(vl,1);
    tpick = 1:min(size(v1,1),vl);
    hold on
    plot(t(ts1(i1)-1+tpick),v1(tpick,i1),':','Color',c{1});
    plot(t(ts2(i2):te2(i2)),v2(1:tl2(i2),i2),':','Color',c{2});
    va(tpick) = va(tpick) + v1(tpick,i1);
    vtpick = 1:tl2(i2);
    tpick = ts2(i2)-ts1(i1)+ vtpick;
    va(tpick) = va(tpick) + v2(vtpick,i2);
    plot(t(ts1(i1):te2(i2)),va,':k');

    tpick = tl(1,i1,i2):tl(2,i1,i2);
    vtpick = tpick - ts1(i1)+1;
    plot(t(tpick),va(vtpick) + b{i1,i2}(:,1),':g');
    
    ax3 = subplot(2,2,4);
    hold on
    plot(t(tl(1,i1,i2):tl(2,i1,i2)),b{i1,i2}(:,2),':k');

    cla(ax1);
    cla(ax2);
    cla(ax3);
end
