function [sEPSP,sIPSP,t] = noAdapVcheck(name,pick,model,picformat,draw,ppp,loadData)
    global tstep printDriver dpi
    msgID = 'MATLAB:rankDeficientMatrix';
    warning('off',msgID);
    visible = true;
    if ~visible
        set(0,'DefaultFigureVisible','off');
    end
    FontSize = 16;
    set(0,'DefaultAxesFontSize',FontSize);
    set(0,'DefaultTextFontSize',FontSize-2);
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
    end
    switch model
        case 'HH'
            silico = @RK4_df;
        case 'IF'
            silico = @RK2_IF_df;
        case 'EIF'
            silico = @RK2_EIF_df;
        otherwise
            disp('model not implemented');
            return
    end
    userpath([pwd,'/channels']);
    name = [name,'-noAdap','-',model];
    load(['parameters-',name]);
    if ~exist(name,'dir')
        mkdir([name,'-vtest']);
    end
    dir = [name,'-vtest'];
%     loadData = false;
    testEE = false;
    testII = false;
    testEI = true;
    testIE = true;
    kQplot = false;
     multipleInput = false;
%    multipleInput = true;
    test = false;
%    test = true;
    ignore = 0;
    ignorefE = 0;
    ignorefI = 0;
    rateE = 0.03;
    rateI = 0.03;
    para.fCurrent = 0;
    meshK = false;
%      ppp = true;
    v0 = -0.4:0.1:1.0;
    v0id = find(abs(v0 - 0.0)<1e-14);
    assert(~isempty(v0id));
    nv0 = length(v0);
    npool = 28;
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if ~isempty(poolobj)
        if poolobj.NumWorkers ~= npool 
            delete(poolobj);
            parpool(npool);
        end
    else
        parpool(npool);
    end
    assert(nv0>1);
    idv0 = [v0id,1:v0id-1,v0id+1:nv0];
    fE = (0.125:0.125:0.5) * 1e-6;
    %fE = (0.025:0.025:0.125) * 1e-6;
    fE = fE(1:4);
    fI = (0.5:0.5:2.0) * 1e-6;
%     fE = (0.5:0.5:1.0) * 1e-6;
%     fI = 0.5*1e-6;
    nE = length(fE);
    nI = length(fI);
    tstep = 0.1;
    %singleStored = false;
    singleStored = true;
    dur = 150;
    dtRange = [0:5:40,50,70,110,150];
    ndtplot = 3;
    tp0 = round(22/tstep);
    t = 0:tstep:dur;
    nt = round(dur/tstep)+1;
    ndt = length(dtRange);
    % 1 2 3
    % E I EI
    i = pick;
    diri = [dir,'/',num2str(i)];
    if ~exist(diri,'dir')
        mkdir(diri);
    end
    vRange = para.vRest(i) + (para.vT(i)-para.vRest(i))*v0;
    if ~loadData
        delete([diri,'/*']);
        if ~singleStored
        [h,sEPSP,sIPSP,vleakage] = sPSP_check(silico,vRange,fI',fE',para,bool,name,dur,i,v0id,dtRange);
            save(['single-',name,'-',num2str(i),'-vTest.mat'],'h','sEPSP','sIPSP','vleakage');
        else
            load(['single-',name,'-',num2str(i),'-vTest.mat']);
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
        if draw
            fname = ['sPSP','-',name,'-',num2str(i)];
            printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
        end
        copyfile('noAdapVcheck.m',[dir,'/noAdapVcheck-',name,'-',num2str(i),'th.m']);
        %
        if draw
            h = figure;
            nsEPSP = squeeze(sEPSP(:,:,1,1,:))./repmat(max(abs(squeeze(sEPSP(:,:,1,1,:)))),[nt,1]);
            nsIPSP = squeeze(sIPSP(:,:,1,1,:))./repmat(max(abs(squeeze(sIPSP(:,:,1,1,:)))),[nt,1]);
            subplot(1,2,1)
            hold on
            c1 = linspace(0,5/6,nv0);
            c2 = linspace(0.3,1,nE)';
            c3 = 0.8;
            for iF = 1:nE
                for iv0 = 1:nv0
                    plot(t,squeeze(nsEPSP(:,iF,iv0)),'Color',hsv2rgb([c1(iv0),c2(iF),c3]));
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
                    plot(t,squeeze(nsIPSP(:,iF,iv0)),'Color',hsv2rgb([c1(iv0),c2(iF),c3]));
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

        kVEE  = zeros(nt,ndt,nv0,nv0);
        kQVEE = zeros(nt,ndt,nv0,nE,nv0);
        kVII  = zeros(nt,ndt,nv0,nv0);
        kQVII = zeros(nt,ndt,nv0,nI,nv0);
        kVEI  = zeros(nt,ndt,nv0,nv0);
        kQVEI = zeros(nt,ndt,nv0,nE,nv0);
        kVIE  = zeros(nt,ndt,nv0,nv0);
        kQVIE = zeros(nt,ndt,nv0,nI,nv0);

        v0case = [1,round(nv0/2),nv0,v0id];

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
        if draw && kQplot
            hQEI = figure;
            hQIE = figure;
            hQEE = figure;
            hQII = figure;
        end
        incVpercent = 10;
        iVpercent = 1;
        for jv0 = 1:nv0
            iv0 = idv0(jv0);
            vPercent = idt/ndt*100;
            if floor(vPercent/incVpercent) > iVpercent
                disp(['--',num2str(round(vPercent)),'%dt']);
                iVpercent = floor(vPercent/incVpercent);
            end
            % redundant definition ignore
%             pv0 = zeros(nE*nI,nt);
%             vadd0 = zeros(nE*nI,nt);
            % % %
            if sum((iv0 - v0case) == 0) && draw
                v0plot = true;
            else
                v0plot = false;
            end
            idtplot = 0;
            incpercent = 10;
            ipercent = 1;
            for idt = 1:ndt
                iidt = round(dtRange(idt)/tstep)+1;
                percent = idt/ndt*100;
                if floor(percent/incpercent) > ipercent
                    disp(['--',num2str(round(percent)),'%dt']);
                    ipercent = floor(percent/incpercent);
                end
                if sum((idt-[1,round(ndt/2),ndt-1])==0)==1 && draw
                    dtplot = true;
                    idtplot = idtplot+1;
                else
                    dtplot = false;
                end
                
                if v0plot && dtplot && ppp
                    pp0 = true;
                else
                    pp0 = false;
                end
                if testEE
                    % EE
                    vv0 = vRange(iv0)*ones(nE*nE,1);
                    para.tE = repmat([0,dtRange(idt)],[nE*nE,1]);
                    para.f_E = [reshape(ones(nE,1)*fE,[nE*nE,1]),repmat(fE',[nE,1])];
                    para.tI = zeros(size(para.tE));
                    para.f_I = zeros(size(para.f_E));
                    [kVEE(:,idt,iv0,:),pv,vadd,vEEDoublet,vEE,iQ,kQ,pQ,h] = doubleCheck(silico,name,para,vv0,bool,tstep,iidt,dur,i,vleakage,sEPSP,sEPSP,vRange,iv0,pp0,v0id,idt,dtRange);
                    kQVEE(:,idt,iv0,:,:) = reshape(kQ,[nt,1,1,nE,nv0]);
                    if v0plot && dtplot
                        hEE0(idtplot) = drawExample(hEE0(idtplot),pv,vadd,pQ,iv0,v0case,vEEDoublet,squeeze(vleakage(:,iv0,idt,:)),t,idt,tp0,kVEE(:,idt,iv0,:),kQ,nE*nE,vRange,para.vRest(i),vEE,'V_{E}V_{E}','EE');
                    end
                    if pp0
                        fname = [num2str(iv0),'v-',num2str(idt),'dt-EE'];
                        printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
                    end
                    if draw && kQplot
                        klim = [-35,35];
                        hQEE = plotkQ(hQEE,t,kQ,nt,iv0,nv0,idt,ndt,klim,iidt,nE);
                    end
                end
                if testII
                    % II
                    vv0 = vRange(iv0)*ones(nI*nI,1);
                    para.tI = repmat([0,dtRange(idt)],[nI*nI,1]);
                    para.f_I = [reshape(ones(nI,1)*fI,[nI*nI,1]),repmat(fI',[nI,1])];
                    para.tE = zeros(size(para.tI));
                    para.f_E = zeros(size(para.f_I));
                    [kVII(:,idt,iv0,:),pv,vadd,vIIDoublet,vII,iQ,kQ,pQ,h] = doubleCheck(silico,name,para,vv0,bool,tstep,iidt,dur,i,vleakage,sIPSP,sIPSP,vRange,iv0,pp0,v0id,idt,dtRange);
                    kQVII(:,idt,iv0,:,:) = reshape(kQ,[nt,1,1,nI,nv0]);
                    if v0plot && dtplot
                        hII0(idtplot) = drawExample(hII0(idtplot),pv,vadd,pQ,iv0,v0case,vIIDoublet,squeeze(vleakage(:,iv0,idt,:)),t,iidt,tp0,kVII(:,idt,iv0,:),kQ,nI*nI,vRange,para.vRest(i),vII,'V_{I}V_{I}','II');
                    end
                    if pp0
                        fname = [num2str(iv0),'v-',num2str(idt),'dt-II'];
                        printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
                    end
                    if draw && kQplot
                        klim = [-35,35];
                        hQII = plotkQ(hQII,t,kQ,nt,iv0,nv0,idt,ndt,klim,iidt,nI);
                    end
                end
                if testEI
                    % EI 
                    vv0 = vRange(iv0)*ones(nE*nI,1);
                    para.f_E = reshape(ones(nI,1)*fE,[nE*nI,1]);
                    para.tE = zeros(nE*nI,1);
                    para.f_I = repmat(fI',[nE,1]);
                    para.tI = dtRange(idt)*ones(nE*nI,1);
                    [kVEI(:,idt,iv0,:),pv,vadd,vEIDoublet,vEI,iQ,kQ,pQ,h] = doubleCheck(silico,name,para,vv0,bool,tstep,iidt,dur,i,vleakage,sEPSP,sIPSP,vRange,iv0,pp0,v0id,idt,dtRange);
                    kQVEI(:,idt,iv0,:,:) = reshape(kQ,[nt,1,1,nE,nv0]);
                    if v0plot && dtplot
                        hEI0(idtplot) = drawExample(hEI0(idtplot),pv,vadd,pQ,iv0,v0case,vEIDoublet,squeeze(vleakage(:,iv0,idt,:)),t,iidt,tp0,kVEI(:,idt,iv0,:),kQ,nE*nI,vRange,para.vRest(i),vEI,'V_{E}V_{I}','EI');
                    end
                    if pp0
                        fname = [num2str(iv0),'v-',num2str(idt),'dt-EI'];
                        printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
                    end
                    if draw && kQplot
                        klim = [-35,35];
                        hQEI = plotkQ(hQEI,t,kQ,nt,iv0,nv0,idt,ndt,klim,iidt,nE);
                    end
                end
                if testIE
                    % IE
                    vv0 = vRange(iv0)*ones(nE*nI,1);
                    para.f_I = reshape(ones(nE,1)*fI,[nE*nI,1]);
                    para.tI = zeros(nE*nI,1);
                    para.f_E = repmat(fE',[nI,1]);
                    para.tE = dtRange(idt)*ones(nE*nI,1);
                    [kVIE(:,idt,iv0,:),pv,vadd,vIEDoublet,vIE,iQ,kQ,pQ,h] = doubleCheck(silico,name,para,vv0,bool,tstep,iidt,dur,i,vleakage,sIPSP,sEPSP,vRange,iv0,pp0,v0id,idt,dtRange);
                    kQVIE(:,idt,iv0,:,:) = reshape(kQ,[nt,1,1,nI,nv0]);
                    if v0plot && dtplot
                        hIE0(idtplot) = drawExample(hIE0(idtplot),pv,vadd,pQ,iv0,v0case,vIEDoublet,squeeze(vleakage(:,iv0,idt,:)),t,iidt,tp0,kVIE(:,idt,iv0,:),kQ,nE*nI,vRange,para.vRest(i),vIE,'V_{I}V_{E}','IE');
                    end
                    if pp0
                        fname = [num2str(iv0),'v-',num2str(idt),'dt-IE'];
                        printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
                    end
                    if draw && kQplot
                        klim = [-35,35];
                        hQIE = plotkQ(hQIE,t,kQ,nt,iv0,nv0,idt,ndt,klim,iidt,nI);
                    end
                end
            end
        end
        assert(idtplot==ndtplot);
        if draw
              %  paper
                height = width/16*16;
                % figure
                mleft = 0.0;    left = mleft*width;
                mbot = 0.0;     bot = mbot*height;    
                mright = 0.0;   fwidth = width-mleft-width*mright;
                mtop = 0.0;     fheight = height-mbot-height*mtop;
                pos = [width, height, left,bot,fwidth,fheight];
                if kQplot
                    if testEE
                        fname = ['kEE','-',name,'-',num2str(i)];
                        printpic(hQEE,diri,fname,picformat,printDriver,dpi,pos);
                    end
                    if testII
                        fname = ['kII','-',name,'-',num2str(i)];
                        printpic(hQII,diri,fname,picformat,printDriver,dpi,pos);
                    end
                    if testEI
                        fname = ['kEI','-',name,'-',num2str(i)];
                        printpic(hQEI,diri,fname,picformat,printDriver,dpi,pos);
                    end
                    if testIE
                        fname = ['kIE','-',name,'-',num2str(i)];
                        printpic(hQIE,diri,fname,picformat,printDriver,dpi,pos);
                    end
                end

            for idt = 1:ndtplot
                if testEE
                    fname = ['EE',num2str(idt),'_example','-',name,'-',num2str(i)];
                    printpic(hEE0(idt),diri,fname,picformat,printDriver,dpi,pos0);
                end
                if testII
                    fname = ['II',num2str(idt),'_example','-',name,'-',num2str(i)];
                    printpic(hII0(idt),diri,fname,picformat,printDriver,dpi,pos0);
                end
                if testEI
                    fname = ['EI',num2str(idt),'_example','-',name,'-',num2str(i)];
                    printpic(hEI0(idt),diri,fname,picformat,printDriver,dpi,pos0);
                end
                if testIE
                    fname = ['IE',num2str(idt),'_example','-',name,'-',num2str(i)];
                    printpic(hIE0(idt),diri,fname,picformat,printDriver,dpi,pos0);
                end
            end
            %if meshK
            %    [X,Y] = meshgrid(0:dtstep:dur,t);
            %    hEI = figure;
            %    subplot(2,1,1)
            %    hold on
            %    c = zeros(nv0,3);
            %    c(:,2) = linspace(0.2,1.0,nv0);
            %    c(:,3) = 1.0;
            %    c = hsv2rgb(c);
            %    
            %    for iv0 = 1:nv0
            %        mesh(X,Y,kEI(:,:,iv0),'EdgeColor',c(iv0,:));
            %    end
            %    xlabel('time (ms)');
            %    ylabel('I(t)-E(t) (ms)');
            %    zlabel('k');
            %    zlim(klim);
            %    subplot(2,1,2)
            %    fname = ['kEImesh','-',name,'-',num2str(i)];
            %    printpic(hEI,diri,fname,picformat,printDriver,dpi,pos0);
            %end
        end
        save([name,'-',num2str(i),'th-vTest'],'kQVEE','kVEE','kQVII','kVII','kQVEI','kVEI','kQVIE','kVIE','tp0','vleakage','sEPSP','sIPSP','dur','vRange','fE','fI','nE','nI','i','dtRange');
    else
        if multipleInput
            load([name,'-',num2str(i),'th-vTest']);
        end
    end
    if multipleInput
        % uniform K
            %>2 input
            % poisson inputs
        %  seed = 237587; % spike
%         poinumber = true;
        poinumber = false;
        poiend = false;
        ndt = ndt-ignore;
        dur0 = dur;
        durpsp = dur;
        nt0 = nt;
        seed =89749;
%         seed = 122435;
        rng(seed);
        l0  = round(durpsp/tstep)+1;
        ldt = round(dtstep/tstep);

        dur = 1000; %ms
        xE = rand(round(rateE*dur*2),1);
        xI = rand(round(rateI*dur*2),1);
        tE = zeros(round(rateE*dur*2),1);
        tI = zeros(round(rateI*dur*2),1);
        pfE = randi(nE-ignorefE,round(rateE*dur*2),1);
        pfI = randi(nI-ignorefI,round(rateI*dur*2),1);
        t = 0;
        it = 0;
        while t < dur
            it = it + 1;
            tE(it) = t;
            t = t - log(xE(it))/rateE;
        end
        para.f_E = fE(pfE(1:it));
        nfE = it;
        para.tE = round(tE(1:it)'/tstep)*tstep;
        t = 0;
        it = 0;
        while t < dur
            it = it +1;
            tI(it) = t;
            t = t - log(xI(it))/rateI;
        end
        para.f_I = fI(pfI(1:it));
        nfI = it;
        para.tI = round(tI(1:it)'/tstep)*tstep;
        fname = ['MultipleInputs',num2str(dur),'-E',num2str(rateE*1e3),'-I',num2str(rateI*1e3),'-S',num2str(seed),'-',name,'-',num2str(i),'th'];
        if test && exist([fname,'.mat'],'file')
            load([fname,'.mat'],'tmp');
        else
            tic;
            para.vtime = -1;
            switch model
                case 'HH'
                    tmp = silico(name,para.vRest(i),para,bool,tstep,dur,i,false);
                case 'EIF'
                    tmp = silico(name,para.vRest(i),para,bool,tstep,dur,i,false,true);
                case 'IF'
                    tmp = silico(name,para.vRest(i),para,bool,tstep,dur,i,false,false);
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

        vE = zeros(nt0,nfE);
        vE0 = zeros(nt0,nfE);
        vE1 = zeros(nt0,nfE);
        vE2 = zeros(nt0,nfE);
        tEs = zeros(nfE,1);
        tEe = zeros(nfE,1);
        tEl = zeros(nfE,1);
        tEeadd = zeros(nfE,1);
        tEladd = zeros(nfE,1);
        vE0pad = zeros(nt,nfE);
        vI = zeros(nt0,nfI);
        vI0 = zeros(nt0,nfI);
        vI1 = zeros(nt0,nfI);
        vI2 = zeros(nt0,nfI);
        tIs = zeros(nfI,1);
        tIe = zeros(nfI,1);
        tIl = zeros(nfI,1);
        tIeadd = zeros(nfI,1);
        tIladd = zeros(nfI,1);
        vI0pad = zeros(nt,nfI);
        vpred0 = zeros(nt,1);
        vpred1 = zeros(nt,1);
        vpred2 = zeros(nt,1);
        vadd = zeros(nt,1);
        disp('linear add:');
        tic; 
        for iE = 1:nfE
            vE0(:,iE) = sEPSP(:,pfE(iE),v0id);
            tEs(iE) = round(para.tE(iE)/tstep)+1;
            tEe(iE) = min(tEs(iE)-1+l0,nt);
            tEl(iE) = tEe(iE)-tEs(iE)+1;
            tEeadd(iE) = min(tEs(iE)-1+nt0,nt);
            tEladd(iE) = tEeadd(iE)-tEs(iE)+1;
            vE0pad(:,iE) = [zeros(tEs(iE)-1,1);vE0(1:tEladd(iE),iE);zeros(nt-tEeadd(iE),1)];
        end
        for iI = 1:nfI
            vI0(:,iI) = sIPSP(:,pfI(iI),v0id);
            tIs(iI) = round(para.tI(iI)/tstep)+1;
            tIe(iI) = min(tIs(iI)-1+l0,nt);
            tIl(iI) = tIe(iI)-tIs(iI)+1;
            tIeadd(iI) = min(tIs(iI)-1+nt0,nt);
            tIladd(iI) = tIeadd(iI)-tIs(iI)+1;
            vI0pad(:,iI) = [zeros(tIs(iI)-1,1);vI0(1:tIladd(iI),iI);zeros(nt-tIeadd(iI),1)];
        end
        vadd0 = sum(vE0pad,2) + sum(vI0pad,2)+ para.vRest(i);
        toc;
        tic;
        vadd = vadd + para.vRest(i);
        iI = 1;
        finishI = false;
        for iE = 1:nfE
            while para.tI(iI) < para.tE(iE)
                vI(:,iI) = interpPSP(sIPSP(:,pfI(iI),:),vadd(tIs(iI)),vRange);
                tpick = tIs(iI):tIeadd(iI);
                vadd(tpick) = vadd(tpick) + vI(1:tIladd(iI),iI);
                if iI == nfI
                    finishI = true;
                    break;
                else
                    iI = iI + 1;
                end
            end
            vE(:,iE) = interpPSP(sEPSP(:,pfE(iE),:),vadd(tEs(iE)),vRange);
            tpick = tEs(iE):tEeadd(iE);
            vadd(tpick) = vadd(tpick) + vE(1:tEladd(iE),iE);
        end
        if ~finishI
            iI0 = iI;
            for iI = iI0:nfI
                vI(:,iI) = interpPSP(sIPSP(:,pfI(iI),:),vadd(tIs(iI)),vRange);
                tpick = tIs(iI):tIeadd(iI);
                vadd(tpick) = vadd(tpick) + vI(1:tIladd(iI),iI);
            end
        end
            %vpred = interpPSPE + interpPSPI;
        toc;
        disp('bilinear correction:');
        tic;
        %vpred0 = vadd0;
        %vpred1 = vadd;
        %vpred2 = vadd;
        bEI0 = cell(nfE,nfI);
        bEI1 = cell(nfE,nfI);
        bEI2 = cell(nfE,nfI);
        tEI = zeros(2,nfE,nfI);  
        bEE0 = cell(nfE,nfE);
        bEE1 = cell(nfE,nfE);
        bEE2 = cell(nfE,nfE);
        tEE = zeros(2,nfE,nfI);
        bIE0 = cell(nfI,nfE);
        bIE1 = cell(nfI,nfE);
        bIE2 = cell(nfI,nfE);
        tIE = zeros(2,nfI,nfE);
        bII0 = cell(nfI,nfI);
        bII1 = cell(nfI,nfI);
        bII2 = cell(nfI,nfI);
        tII = zeros(2,nfI,nfI);

        vpred0 = vpred0 + para.vRest(i);
        vpred1 = vpred1 + para.vRest(i);
        vpred2 = vpred2 + para.vRest(i);

        finishI = false;
        iI = 1;
        %debug = true;
        debug = false;
        if debug
            dbstop if error
            dbstop at 1284
            hDebug = figure;
            xDebug = [415,480];
            hv0 = gobjects(1);
            hv1 = gobjects(1);
            hv2 = gobjects(1);
            subplot(2,1,1);
            hold on
            plot(t,v,'k');
            plot(t,vadd0,':k');
            plot(t,vadd,':m');
        end
          rule = @(iiE,iiI,iE,iI) iiI == 10 && iE == 26;
        rule = @(iiE,iiI,iE,iI) iI == 10;
        for iE = 1:nfE
            while para.tI(iI) < para.tE(iE) && ~finishI
                %vI1(:,iI) = vI(:,iI);
                %vI2(:,iI) = vI(:,iI);
                vI1(:,iI) = interpPSP(sIPSP(:,pfI(iI),:),vpred1(tIs(iI)),vRange);
                vI2(:,iI) = interpPSP(sIPSP(:,pfI(iI),:),vpred2(tIs(iI)),vRange);
                tpick = tIs(iI):tIeadd(iI);
                vpred0(tpick) = vpred0(tpick) + vI0(1:tIladd(iI),iI);
                vpred1(tpick) = vpred1(tpick) + vI1(1:tIladd(iI),iI);
                vpred2(tpick) = vpred2(tpick) + vI2(1:tIladd(iI),iI);
                if testEI
                    for iiE = 1: iE-1
                        if para.tI(iI) - para.tE(iiE) < durpsp
                            if debug && rule(iiE,0,0,iI)
                                hDebug = debugPreplot(hDebug,t,vpred0,vpred1,vpred2);
                            end
                            [vpred0, vpred1, vpred2, bEI0{iiE,iI},bEI1{iiE,iI},bEI2{iiE,iI},tEI(:,iiE,iI)] = pred(vpred0,vpred1,vpred2,vadd,vadd,...
                                                                                                                  iiE,iI,vE0,vE1,vE2,vI0,vI1,vI2,...
                                                                                                                  tEs,tEe,tEl,tIs,kEI0,kEI,kQEI,pfE,...
                                                                                                                  vRange,nv0,ndt,ldt);
                            if debug && rule(iiE,0,0,iI)
                                hDebug = plotDebug(hDebug,t,v,vpred0,vpred1,vpred2,iiE,iI,...
                                                   bEI0,bEI1,bEI2,tEI,pfE,pfI,...
                                                   tEs,tEe,tEl,tIs,tIe,tIl,'EI',xDebug,...
                                                   vE0,vE1,vE2,vI0,vI1,vI2);
                                rule = @(iiE,iiI,iE,iI) true;
                            end
                        end
                    end
                end
                if testII
                    for iiI = 1: iI-1
                        if para.tI(iI) - para.tI(iiI) < durpsp
                            if debug && rule(0,iiI,0,iI)
                                hDebug = debugPreplot(hDebug,t,vpred0,vpred1,vpred2);
                            end
                            [vpred0, vpred1, vpred2, bII0{iiI,iI},bII1{iiI,iI},bII2{iiI,iI},tII(:,iiI,iI)] = pred(vpred0,vpred1,vpred2,vadd,vadd,...
                                                                                                                  iiI,iI,vI0,vI1,vI2,vI0,vI1,vI2,...
                                                                                                                  tIs,tIe,tIl,tIs,kII0,kII,kQII,pfI,...
                                                                                                                  vRange,nv0,ndt,ldt);
                            if debug && rule(0,iiI,0,iI)
                                hDebug = plotDebug(hDebug,t,v,vpred0,vpred1,vpred2,iiI,iI,...
                                                   bII0,bII1,bII2,tII,pfI,pfI,...
                                                   tIs,tIe,tIl,tIs,tIe,tIl,'II',xDebug,...
                                                   vI0,vI1,vI2,vI0,vI1,vI2);
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
            %vE1(:,iE) = vE(:,iE);
            %vE2(:,iE) = vE(:,iE);

            vE1(:,iE) = interpPSP(sEPSP(:,pfE(iE),:),vpred1(tEs(iE)),vRange);
            vE2(:,iE) = interpPSP(sEPSP(:,pfE(iE),:),vpred2(tEs(iE)),vRange);
            tpick = tEs(iE):tEeadd(iE);
            vpred0(tpick) = vpred0(tpick) + vE0(1:tEladd(iE),iE);
            vpred1(tpick) = vpred1(tpick) + vE1(1:tEladd(iE),iE);
            vpred2(tpick) = vpred2(tpick) + vE2(1:tEladd(iE),iE);
            if testIE
                for iiI = 1:iI-1
                    if para.tE(iE) - para.tI(iiI) < durpsp
                        if debug && rule(0,iiI,iE,0)
                            hDebug = debugPreplot(hDebug,t,vpred0,vpred1,vpred2);
                        end
                        [vpred0, vpred1, vpred2, bIE0{iiI,iE},bIE1{iiI,iE},bIE2{iiI,iE},tIE(:,iiI,iE)] = pred(vpred0,vpred1,vpred2,vadd,vadd,...
                                                                                                              iiI,iE,vI0,vI1,vI2,vE0,vE1,vE2,...
                                                                                                              tIs,tIe,tIl,tEs,kIE0,kIE,kQIE,pfI,...
                                                                                                              vRange,nv0,ndt,ldt);
                        if debug && rule(0,iiI,iE,0)
                            hDebug = plotDebug(hDebug,t,v,vpred0,vpred1,vpred2,iiI,iE,...
                                               bIE0,bIE1,bIE2,tIE,pfI,pfE,...
                                               tIs,tIe,tIl,tEs,tEe,tEl,'IE',xDebug,...
                                               vI0,vI1,vI2,vE0,vE1,vE2);
                            rule = @(iiE,iiI,iE,iI) true;
                        end
                    end
                end
            end
            if testEE
                for iiE = 1:iE-1
                    if para.tE(iE) - para.tE(iiE) < durpsp
                        if debug && rule(iiE,0,iE,0)
                            hDebug = debugPreplot(hDebug,t,vpred0,vpred1,vpred2);
                        end
                        [vpred0, vpred1, vpred2, bEE0{iiE,iE},bEE1{iiE,iE},bEE2{iiE,iE},tEE(:,iiE,iE)] = pred(vpred0,vpred1,vpred2,vadd,vadd,...
                                                                                                              iiE,iE,vE0,vE1,vE2,vE0,vE1,vE2,...
                                                                                                              tEs,tEe,tEl,tEs,kEE0,kEE,kQEE,pfE,...
                                                                                                              vRange,nv0,ndt,ldt);
                        if debug && rule(iiE,0,iE,0)
                            hDebug = plotDebug(hDebug,t,v,vpred0,vpred1,vpred2,iiE,iE,...
                                               bEE0,bEE1,bEE2,tEE,pfE,pfE,...
                                               tEs,tEe,tEl,tEs,tEe,tEl,'EE',xDebug,...
                                               vE0,vE1,vE2,vE0,vE1,vE2);
                            rule = @(iiE,iiI,iE,iI) true;
                        end
                    end
                end
            end
        end
        if ~finishI
            iI0 = iI;
            for iI = iI0:nfI
                %vI1(:,iI) = vI(:,iI);
                %vI2(:,iI) = vI(:,iI);

                vI1(:,iI) = interpPSP(sIPSP(:,pfI(iI),:),vpred1(tIs(iI)),vRange);
                vI2(:,iI) = interpPSP(sIPSP(:,pfI(iI),:),vpred2(tIs(iI)),vRange);
                tpick = tIs(iI):tIeadd(iI);
                vpred0(tpick) = vpred0(tpick) + vI0(1:tIladd(iI),iI);
                vpred1(tpick) = vpred1(tpick) + vI1(1:tIladd(iI),iI);
                vpred2(tpick) = vpred2(tpick) + vI2(1:tIladd(iI),iI);
                if testEI
                    for iiE = 1: nfE
                        if para.tI(iI) - para.tE(iiE) < durpsp
                            if debug && rule(iiE,0,0,iI)
                                hDebug = debugPreplot(hDebug,t,vpred0,vpred1,vpred2);
                            end
                            [vpred0, vpred1, vpred2, bEI0{iiE,iI},bEI1{iiE,iI},bEI2{iiE,iI},tEI(:,iiE,iI)] = pred(vpred0,vpred1,vpred2,vadd,vadd,...
                                                                                                                  iiE,iI,vE0,vE1,vE2,vI0,vI1,vI2,...
                                                                                                                  tEs,tEe,tEl,tIs,kEI0,kEI,kQEI,pfE,...
                                                                                                                  vRange,nv0,ndt,ldt);
                            if debug && rule(iiE,0,0,iI)
                                hDebug = plotDebug(hDebug,t,v,vpred0,vpred1,vpred2,iiE,iI,...
                                                   bEI0,bEI1,bEI2,tEI,pfE,pfI,...
                                                   tEe,tEl,tIs,tIe,tIl,'EI',xDebug,...
                                                   vE1,vE2,vI0,vI1,vI2);
                                rule = @(iiE,iiI,iE,iI) true;
                            end
                        end
                    end
                end
                if testII
                    for iiI = 1: iI-1
                        if para.tI(iI) - para.tI(iiI) < durpsp
                            if debug && rule(0,iiI,0,iI)
                                hDebug = debugPreplot(hDebug,t,vpred0,vpred1,vpred2);
                            end
                            [vpred0, vpred1, vpred2, bII0{iiI,iI},bII1{iiI,iI},bII2{iiI,iI},tII(:,iiI,iI)] = pred(vpred0,vpred1,vpred2,vadd,vadd,...
                                                                                                                  iiI,iI,vI0,vI1,vI2,vI0,vI1,vI2,...
                                                                                                                  tIs,tIe,tIl,tIs,kII0,kII,kQII,pfI,...
                                                                                                                  vRange,nv0,ndt,ldt);
                            if debug && rule(0,iiI,0,iI)
                                hDebug = plotDebug(hDebug,t,v,vpred0,vpred1,vpred2,iiI,iI,...
                                                   bII0,bII1,bII2,tII,pfI,pfI,...
                                                   Ie,tIl,tIs,tIe,tIl,'II',xDebug,...
                                                   I1,vI2,vI0,vI1,vI2);
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
        p = plot(t,v);
        plot(t,vadd0,':','Color',p.Color);
        plot(t,vpred0,'--','Color',p.Color);
        plot(t,vadd,':m');
        plot(t,vpred1,'--m');
        plot(t,vpred2,'--g');
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
            plot(ones(2,1)*para.tE(selectE)+durpsp,[v(tEs(selectE)+l0-1)';y(2)*ones(1,sfE)],':r');
            plot(ones(2,1)*para.tI(selectI)+durpsp,[v(tIs(selectI)+l0-1)';y(2)*ones(1,sfI)],':b');
        end
        if poinumber
            for iE = 1:nfE
                if para.tE(iE) > xl(1) && para.tE(iE) < xl(2)
                    text(para.tE(iE),vTargetE(iE)-sl*0.1,num2str(iE),'Color','r','FontSize',FS);
                    text(para.tE(iE),vTargetE(iE)-sl*0.2,num2str(pfE(iE)),'Color','r','FontSize',FS);
                end
                if poiend
                    if para.tE(iE)+durpsp > xl(1) && para.tE(iE)+durpsp < xl(2)
                        text(para.tE(iE)+durpsp,v(tEs(iE)+l0-1)+sl*0.1,num2str(iE),'Color','r','FontSize',FS);
                        text(para.tE(iE)+durpsp,v(tEs(iE)+l0-1)+sl*0.2,num2str(pfE(iE)),'Color','r','FontSize',FS);
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
                        text(para.tI(iI)+durpsp,v(tIs(iI)+l0-1)+sl*0.1,num2str(iI),'Color','b','FontSize',FS);
                        text(para.tI(iI)+durpsp,v(tIs(iI)+l0-1)+sl*0.2,num2str(pfI(iI)),'Color','b','FontSize',FS);
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
        errorbar(0.95,sign(sum(vadd0-v>0)*2-1) *mean(abs(vadd0-v)), std(abs(vadd0-v)),'*b');
        errorbar(1.05,sign(sum(vadd -v>0)*2-1) *mean(abs(vadd-v)),  std(abs(vadd-v)),'*m');
        e1=errorbar(1.9,sign(sum(vpred0-v>0)*2-1)*mean(abs(vpred0-v)),std(abs(vpred0-v)),'*b');
        e2=errorbar(2.0,sign(sum(vpred1-v>0)*2-1)*mean(abs(vpred1-v)),std(abs(vpred1-v)),'*m');
        e3=errorbar(2.1,sign(sum(vpred2-v>0)*2-1)*mean(abs(vpred2-v)),std(abs(vpred2-v)),'*g');
        set(ax3,'XTick',[1,2]);
        set(ax3,'XTickLabel',{'linear';'bilinear'});
        legend([e1,e2,e3],{'K_{0}','K','K_{Q}'});
        ylabel('error mV');
        xlim([0.5,2.5]);
        ax4 = subplot(2,2,4);
%         xl = [t(min(tEs(round(nfE/2)),tIs(round(nfI/2)))),t(nt)];
        ix = 23;
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
                text(para.tE(iE)+durpsp,v(tEs(iE)+l0-1)+sl*0.1,num2str(iE),'Color','r','FontSize',FS);
                text(para.tE(iE)+durpsp,v(tEs(iE)+l0-1)+sl*0.2,num2str(pfE(iE)),'Color','r','FontSize',FS);
            end
        end
        for iI = 1:nfI
            if para.tI(iI) > xl(1) && para.tI(iI) < xl(2)
                text(para.tI(iI),vTargetI(iI)-sl*0.1,num2str(iI),'Color','b','FontSize',FS);
                text(para.tI(iI),vTargetI(iI)-sl*0.2,num2str(pfI(iI)),'Color','b','FontSize',FS);
            end
            if para.tI(iI)+durpsp > xl(1) && para.tI(iI)+durpsp < xl(2)
                text(para.tI(iI)+durpsp,v(tIs(iI)+l0-1)+sl*0.1,num2str(iI),'Color','b','FontSize',FS);
                text(para.tI(iI)+durpsp,v(tIs(iI)+l0-1)+sl*0.2,num2str(pfI(iI)),'Color','b','FontSize',FS);
            end
        end
        if ~poiend
            selectE = (para.tE + durpsp < xl(2) ) & (para.tE + durpsp > xl(1));
            sfE = sum(selectE);
            selectI = (para.tI + durpsp < xl(2) ) & (para.tI + durpsp > xl(1));
            sfI = sum(selectI);
            plot(ones(2,1)*para.tE(selectE)+durpsp,[v(tEs(selectE)+l0-1)';y(2)*ones(1,sfE)],':r');
            plot(ones(2,1)*para.tI(selectI)+durpsp,[v(tIs(selectI)+l0-1)';y(2)*ones(1,sfI)],':b');
        end
        xlabel('ms');
        ylabel('mV');
        xlim(xl);
        
        ylim([y2(1),y2(2)]);
        printpic(hM,diri,fname,picformat,printDriver,dpi,pos0);
    end
end
function [h,EPSP,IPSP,vleakage] = sPSP_check(silico,vRange,fIRange,fERange,para,bool,name,dur,i,v0id,dtRange)
    global tstep
    ndt = length(dtRange);
    nv0 = length(vRange);
    nE = length(fERange);
    nI = length(fIRange);
    nt = round(dur/tstep)+1;
    EPSP = zeros(nt,nE,nv0,ndt,nv0);
    IPSP = zeros(nt,nI,nv0,ndt,nv0);
    vleakage = zeros(nt,nv0,ndt,nv0);

    v0 = vRange';
    para.f_E = zeros(nv0,1);
    para.tE = para.f_E;
    para.f_I = zeros(nv0,1);
    para.tI = para.f_I;
    for iv0 = 1:nv0
        para.newv = vRange(iv0);
        parfor idt = 1:ndt
            param = para;
            param.vtime = dtRange(idt);
            leakage = silico(name,v0,param,bool,tstep,dur,i,false);
            vleakage(:,:,idt,iv0) = squeeze(leakage(1,:,:));
        end
    end
    disp('leakage complete');
    for iv0 = 1:nv0
        disp([num2str(iv0),'th v complete']);
        for iiv0 = 1:nv0
            para.newv = vRange(iiv0);

            v0 = vRange(iv0)*ones(nE,1);
            para.f_E = fERange;
            para.tE = zeros(length(fERange),1);
            para.f_I = para.tE;
            para.tI = para.tE;

            parfor idt = 1:ndt
                param = para;
                param.vtime = dtRange(idt);
                tmpE = silico(name,v0,param,bool,tstep,dur,i,false);
                EPSP(:,:,iv0,idt,iiv0) = squeeze(tmpE(1,:,:)) - repmat(vleakage(:,iv0,idt,iiv0),[1,nE]);
            end

            v0 = vRange(iv0)*ones(nI,1);
            para.f_I = fIRange;
            para.tI = zeros(length(fIRange),1);
            para.f_E = para.tI;
            para.tE = para.tI;
             
            parfor idt = 1:ndt
                param = para;
                param.vtime = dtRange(idt);
                tmpI = silico(name,v0,param,bool,tstep,dur,i,false);
                IPSP(:,:,iv0,idt,iiv0) = squeeze(tmpI(1,:,:)) - repmat(vleakage(:,iv0,idt,iiv0),[1,nI]);
            end
        end
    end
    idt = 1;
    h = figure;
    s = [0.1,1.0];
    v = 1.0;
    t = linspace(0,dur,nt);
    subplot(2,2,2);
    plotSinglePSP(0,s,v,fERange,t,squeeze(EPSP(:,:,2,idt,:)),nE,nv0,v0id);
    xlabel('f');
    ylabel('t');
    zlabel('EPSP (mV)');
    subplot(2,2,4);
    plotSinglePSP(2/3,s,v,fIRange,t,squeeze(IPSP(:,:,2,idt,:)),nI,nv0,v0id);
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
        plot(t,squeeze(vleakage(:,1,idt,:)),'Color',c(iv0,:));
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
function [kV,pV,vaddV,vDoubletV,EI,iQV,kQV,pQV,h] = doubleCheck(silico,name,para,v0,bool,tstep,iidt,dur,i,vleakage,sPSP1,sPSP2,vRange,iv0,pp,v0id,idt,dtRange)
     n1 = size(sPSP1,2);
     n2 = size(sPSP2,2);
     nt = round(dur/tstep)+1;
     nv0 = length(vRange);

     kV = zeros(nt,nv0);
     vDoubletV = zeros(nt,n1*n2,nv0);
     vaddV = zeros(n1*n2,nt,nv0);
     v1V = zeros(n1*n2,nt,nv0);
     v2V = zeros(n1*n2,nt,nv0);
     interpPSP2V = zeros(nt,n1*n2,nv0);
     vt1V = zeros(nt,n1*n2,nv0);
     vt2V = zeros(nt,n1*n2,nv0);
     pV = zeros(n1*n2,nt,nv0);
     pQV = pV;
     pQV0 = pV;
     iQV = zeros(n1,nv0);
     kQV = zeros(nt,n1,nv0);
     para.vtime = (iidt-1)*tstep;
     parfor iiv0 = 1:nv0
        param = para;
        param.newv = vRange(iiv0);
        tmp = silico(name,v0,param,bool,tstep,dur,i,false);
        vDoubletV(:,:,iiv0) = squeeze(tmp(1,:,:));
        vDoubletV(:,:,iiv0) = vDoubletV(:,:,iiv0) - repmat(vleakage(:,iv0,idt,iiv0),[1,n1*n2]);
     end
     for iiv0 = 1:nv0
        vTarget = vRange(iiv0);
        for iF = 1:n1
            range = (iF-1)*n2+(1:n2);
            v1V(range,:,iiv0) = (sPSP1(:,iF,iv0,idt,iiv0)*ones(1,n2))';
            [interpPSP2V(:,range,iiv0), vt1V(:,range,iiv0),vt2V(:,range,iiv0), ~] = interpPSP(sPSP2(:,:,1,1,:),vTarget,vRange);
            [~, tmp] = max(abs(interpPSP2V(1:nt-iidt+1,range,iiv0)));
            iQV(iF,iiv0) = iidt -1 + round(mean(tmp));
            %disp(std(tmp));
            v2V(range,:,iiv0) = [zeros(iidt-1,n2);interpPSP2V(1:nt-iidt+1,range,iiv0)]';
            vaddV(range,:,iiv0) = v1V(range,:,iiv0) + v2V(range,:,iiv0);
            parfor it = 1:nt
               [kQV(it,iF,iiv0), pQV(range,it,iiv0)] = p_fit110k0(v1V(range,it,iiv0),v2V(range,it,iiv0),vDoubletV(it,range,iiv0)');
            end
            pQV0(range,:,iiv0) = vaddV(range,:,iiv0) + kQV(iQV(iF,iiv0),iF,iiv0)*v1V(range,:,iiv0).*v2V(range,:,iiv0);
        end
        parfor it = 1:nt
            [kV(it,iiv0),pV(:,it,iiv0)] = p_fit110k0(v1V(:,it,iiv0),v2V(:,it,iiv0),vDoubletV(it,:,iiv0)');
        end
     end

     EI = v1V.*v2V;
     h = 0;
     if pp 
        h = figure;
        t = 0:tstep:dur;
        subplot(3,n2,1)
        text(-0.2,1.2,[num2str(iv0),'th v0, linear sum in dot'],'Units','Normalized');
        subplot(3,n2,1*n2+1)
        text(-0.2,1.2,'Zoom, predict in dash','Units','Normalized');
        subplot(3,n2,2*n2+1)
        text(-0.2,1.2,'Interpolation','Units','Normalized');
        ii = 1;
        iiv0 = round(nv0/2);
        range = (ii-1)*n2+(1:n2);
        for i2 =1:n2
            subplot(3,n2,i2);
            title(['IPSP ',num2str(i2)]);
            hold on
            plot(t,v1V(range(1),:,iiv0),':r');
            plot(t,v2V(range(i2),:,iiv0),':b');
            plot(t,vaddV(range(i2),:,iiv0),':k');
            plot(t,vDoubletV(:,range(i2),iiv0),'k');
            subplot(3,n2,n2+i2);
            hold on
            pick = iidt:nt;
            plot(t(pick),v1V(range(1),pick,iiv0),':r');
            plot(t(pick),v2V(range(i2),pick,iiv0),':b');
            plot(t(pick),vaddV(range(i2),pick,iiv0),':k');
            plot(t(pick),pV(range(i2),pick,iiv0),'--k');
            plot(t(pick),pQV0(range(i2),pick,iiv0),':m');
            plot(t(pick),pQV(range(i2),pick,iiv0),'--m');
            [ax,h1,h2] = plotyy(t(pick),vDoubletV(pick,range(i2),iiv0),t(pick),kV(pick,iiv0));
            h1.Color = 'k';
            h2.LineStyle = ':';
            h2.Color = 'g';
            hold(ax(2),'on');
            plot(ax(2),t(pick),kQV(pick,ii,iiv0),'--g');
            plot(ax(2),t(iQV(ii,iiv0)),kQV(iQV(ii,iiv0),ii,iiv0),'*g');
            ylim(ax(2),[-100,100]);
            subplot(3,n2,2*n2+i2);
            hold on
            pick2 = 1:nt-iidt+1;
            plot(t(pick),vt1V(pick2,range(i2),iiv0),':r');
            plot(t(pick),vt2V(pick2,range(i2),iiv0),':b');
            plot(t(pick),interpPSP2V(pick2,range(i2),iiv0),':k');
        end
     end
end
function h = drawExample(h,pv0,vadd0,pQ0,iv0,v0case,vDoublet0,vleakage0,t,iidt,tp0,k0,kQ0,n,vRange,vRest,v120,xl,titl)
    nCase = length(v0case);
    jv0 = find((iv0 - v0case) ==0, 1,'first');
    nv0 = length(vRange);
    inv0 = 0;
    for iiv0 = [1,round(nv0/2),nv0]
        vleakage = repmat(vleakage0(:,iiv0),[1,n]);
        k = k0(:,iiv0);
        v12 = v120(:,:,iiv0);
        kQ = kQ0(:,:,iiv0);
        vDoublet = vDoublet0(:,:,iiv0) + vleakage - vRest;
        pv = pv0(:,:,iiv0) + vleakage' - vRest;
        pQ = pQ0(:,:,iiv0) + vleakage' - vRest;
        vadd = vadd0(:,:,iiv0) + vleakage' - vRest;
        nt = length(t);
        figure(h);
        tp = min(iidt + tp0,nt);
        subplot(5,nCase,jv0);
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
        title([titl,', v0 = ',num2str(vRange(iv0),'%3.1f')]);

        subplot(5,nCase,nCase+jv0);
        hold on
        plot(vDoublet(tp,:)',pv(:,tp),'o');
        plot(vDoublet(tp,:)',vadd(:,tp),'s');
        plot(vDoublet(tp,:)',pQ(:,tp),'d');
        x = xlim;
        y = x;
        plot(x,y,':k');
        xlabel('PSP');
        ylabel('predict');

        subplot(5,nCase,nCase*(2+inv0)+jv0)
        hold on
        p = plot(v12(:,tp),(vDoublet(tp,:)'-vadd(:,tp)),'o');
        xx = sort(v12(:,tp));
        plot(xx,xx*k(tp),':k','LineWidth',1.5);
        plot(xx,(xx*kQ(tp,:)),'--','Color',p.Color);
        ylabel('V_{SC}');
        xlabel(xl);
        title({['k= ',num2str(k(tp),'%.3g')]});
        inv0 = inv0 + 1;
    end
end
function hQ = plotkQ(hQ,t,kQ,nt,iv0,nv0,idt,ndt,klim,iidt,n1)
        figure(hQ)
%         nkQ = kQ./repmat(max(abs(kQ),[],2),[1,nt]);
        subplot(ndt,2,(idt-1)*2+1)
        if idt == 1
            title('hue v0, saturation f, \Delta_{tInput} subplots');
        end
        hold on
        if iv0 == 1
            plot([t(iidt),t(iidt)],klim,':k');
        end 
        for iF = 1:n1
            color = hsv2rgb([0+(iv0-1)/(nv0-1)*5/6,0.2+(iF-1)/(n1-1)*0.8,0.8]);
            plot(t,kQ(iF,:),'Color',color);
        end
        if idt ~= ndt
            set(gca,'XTickLabel','');
        end
        ylim(klim);

        subplot(nv0,2,iv0*2)
        if iv0 == 1
            title('hue dt, saturation f, v0 subplots');
        end
        hold on
        plot([t(iidt),t(iidt)],klim,':k');
        for iF = 1:n1
            color = hsv2rgb([0+(idt-1)/(ndt-1)*5/6,0.3+(iF-1)/(n1-1)*0.7,0.8]);
            plot(t(iidt:nt),kQ(iF,iidt:nt),'Color',color);
        end
        if iv0 ~= nv0
            set(gca,'XTickLabel','');
        end
        ylim(klim);
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
function [vpred0, vpred1, vpred2, b0, b1, b2, t] = pred(vpred0,vpred1,vpred2,...
                                                        v1,v2,i1,i2,v10,v11,v12,v20,v21,v22,...
                                                        t1s,t1e,t1l,t2s,k0,k,kQ,pf,vRange,nv0,ndt,ldt)
    lt = t1e(i1)-t2s(i2);
    assert(lt>0);
    tpick = t2s(i2)+(0:lt);
    t(1) = tpick(1);
    t(2) = tpick(end);
    b0 = zeros(length(tpick),2);
    b1 = b0;
    b2 = b0;
    rdt = (t1l(i1)-1-lt)/ldt;
    fdt = floor(rdt);
    cdt = fdt+1;
    k0tmp = interpK0(rdt,fdt,cdt,k0,ndt,ldt);
    b0(:,1) = k0tmp(1:lt+1).*v10((t1l(i1)-lt):t1l(i1),i1).*v20(1:lt+1,i2);
    b0(:,2) = k0tmp(1:lt+1);
    vpred0(tpick) = vpred0(tpick) + b0(:,1);

    k1tmp = interpK(rdt,fdt,cdt,v1(t1s(i1)),vRange,nv0,k,ndt,ldt);
    b1(:,1) = k1tmp(1:lt+1).*v11((t1l(i1)-lt):t1l(i1),i1).*v21(1:lt+1,i2);
    b1(:,2) = k1tmp(1:lt+1);
    vpred1(tpick) = vpred1(tpick) + b1(:,1);

    k2tmp = interpK(rdt,fdt,cdt,v2(t1s(i1)),vRange,nv0,kQ(:,:,:,pf(i1)),ndt,ldt);
    b2(:,1) = k2tmp(1:lt+1).*v12((t1l(i1)-lt):t1l(i1),i1).*v22(1:lt+1,i2);
    b2(:,2) = k2tmp(1:lt+1);
    vpred2(tpick) = vpred2(tpick) + b2(:,1);
end


function k0tmp = interpK0(rdt,fdt,cdt,k,ndt,ldt)
    k0 = k( (fdt*ldt+1) : ((ndt-1)*ldt+1), fdt+1);
    k0tmp = k0 + ...
            (rdt-fdt) * (k( (cdt*ldt+1):(ndt*ldt+1), cdt+1) - k0);
end
function ktmp = interpK(rdt,fdt,cdt,vTarget,v,nv0,k,ndt,ldt)
    t0 = (fdt*ldt+1) : ((ndt-1)*ldt+1);
    t1 = (cdt*ldt+1) : (ndt*ldt+1);
    fv = find(vTarget-v>=0, 1,'last');
    if isempty(fv)
        k0 = k(t0,fdt+1,1);
        ktmp = k0 + (vTarget-v(1))/(v(1)-v(2)) * (k0-k(t0,fdt+1,2)) + (rdt-fdt) * (k(t1,cdt+1,1)-k0);
    else
        if fv == nv0
            k0 = k(t0,fdt+1,fv);
            ktmp = k0 + (vTarget-v(fv))/(v(fv)-v(fv-1)) * (k0-k(t0,fdt+1,fv-1)) + (rdt-fdt) * (k(t1,cdt+1,fv)-k0);
        else
            k0 = k(t0,fdt+1,fv);
            ktmp = k0 + (vTarget-v(fv))/(v(fv+1)-v(fv))*(k(t0,fdt+1,fv+1)-k0) + (rdt-fdt) * (k(t1,cdt+1,fv)-k0);
        end
    end
%     figure;
%     hold on
%     plot(k0,'b')
%     plot(k(t0,fdt+1,fv+1),'r');
%     plot(k(t0,cdt+1,fv),':g');
%     plot(k(t1,cdt+1,fv),'g');
%     plot(ktmp,'k');
%     title(num2str(rdt-fdt));
%     return
end
function h = debugPreplot(h,t,vpred0,vpred1,vpred2)
    figure(h)
    subplot(2,1,1)
    hold on
    plot(t,vpred0,'-b');
    plot(t,vpred1,'-m');
    plot(t,vpred2,'-g');
end
function h = plotDebug(h,t,v,vpred0,vpred1,vpred2,i1,i2,...
                       b0,b1,b2,tl,f1,f2,...
                       ts1,te1,tl1,ts2,te2,tl2,label,xl,...
                       v10,v11,v12,v20,v21,v22)
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
    plot(t,vpred0,'--b');
    plot(t,vpred1,'--m');
    plot(t,vpred2,'--g');
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
    va0 = zeros(vl,1);
    va1 = zeros(vl,1);
    va2 = zeros(vl,1);
    tpick = 1:min(size(v10,1),vl);
    hold on
    plot(t(ts1(i1)-1+tpick),v10(tpick,i1),':','Color',c{1});
    plot(t(ts2(i2):te2(i2)),v20(1:tl2(i2),i2),':','Color',c{2});
    plot(t(ts1(i1)-1+tpick),v11(tpick,i1),'--','Color',c{1});
    plot(t(ts2(i2):te2(i2)),v21(1:tl2(i2),i2),'--','Color',c{2});
    plot(t(ts1(i1)-1+tpick),v12(tpick,i1),'-','Color',c{1});
    plot(t(ts2(i2):te2(i2)),v22(1:tl2(i2),i2),'-','Color',c{2});
    va0(tpick) = va0(tpick) + v10(tpick,i1);
    va1(tpick) = va1(tpick) + v11(tpick,i1);
    va2(tpick) = va2(tpick) + v12(tpick,i1);
    vtpick = 1:tl2(i2);
    tpick = ts2(i2)-ts1(i1)+ vtpick;
    va0(tpick) = va0(tpick) + v20(vtpick,i2);
    va1(tpick) = va1(tpick) + v21(vtpick,i2);
    va2(tpick) = va2(tpick) + v22(vtpick,i2);
    plot(t(ts1(i1):te2(i2)),va0,':k');
    plot(t(ts1(i1):te2(i2)),va1,'--k');
    plot(t(ts1(i1):te2(i2)),va2,'-k');

    tpick = tl(1,i1,i2):tl(2,i1,i2);
    vtpick = tpick - ts1(i1)+1;
    plot(t(tpick),va0(vtpick) + b0{i1,i2}(:,1),':g');
    plot(t(tpick),va1(vtpick) + b1{i1,i2}(:,1),'--g');
    plot(t(tpick),va2(vtpick) + b2{i1,i2}(:,1),'-g');
    
    ax3 = subplot(2,2,4);
    hold on
    plot(t(tl(1,i1,i2):tl(2,i1,i2)),b0{i1,i2}(:,2),':k');
    plot(t(tl(1,i1,i2):tl(2,i1,i2)),b1{i1,i2}(:,2),'--k');
    plot(t(tl(1,i1,i2):tl(2,i1,i2)),b2{i1,i2}(:,2),'-k');

    cla(ax1);
    cla(ax2);
    cla(ax3);
end
