function [sEPSP,sIPSP,t] = noAdapV_k4(theme,name,pick,model,picformat,draw,ppp,loadData,npool,v0,fE,fI,singleStored, dur, dtRange, tstep, rateE, rateI, mdur,linear0,bilinear0)
    global nvplot ndtplot iv0Case idtCase
    msgID = 'MATLAB:rankDeficientMatrix';
    warning('off',msgID);
    visible = true;
    if ~visible
        set(0,'DefaultFigureVisible','off');
    end
    FontSize = 16;
    set(0,'DefaultAxesFontSize',FontSize);
    set(0,'DefaultTextFontSize',FontSize-2);
    if nargin < 20
        bilinear0 = false;
        if nargin < 19
            linear0 = false;
            if nargin < 18
                rateI = 20;
                if nargin < 17
                    rateE = 40;
                    if nargin < 16
                        tstep = 1.0/32;
                        if nargin < 15
                            dtRange = -12;
                            if nargin < 14
                                dur = 300;
                                if nargin < 13
                                    singleStored = false;
                                    if nargin < 12
                                        fI = linspace(0.5,2.0,4) * 1e-5;
                                        if nargin < 11
                                            fE = linspace(0.125,0.5,4) * 1e-5;
                                            if nargin < 10
                                                v0 = -0.4:0.1:1.2;
                                                if nargin < 9
                                                    npool = 1;
                                                    if nargin < 8
                                                        loadData = true;
                                                        if nargin < 7
                                                            ppp = false;
                                                            if nargin < 6
                                                                draw = false;
                                                                if nargin < 5
                                                                    picformat = '';
                                                                    if nargin < 4
                                                                        model = 'HH';
                                                                        if nargin < 3
                                                                            pick = 1;
                                                                        end
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    explore = false;
    loadExplore = false;
    if dur < 0 || dtRange < 0
        explore = true;
    end
    if dur == 0 && dtRange == 0
        loadExplore = true;
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
    picformat
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
    pname = ['../library/parameters-',name,'-noAdap-',model];
    disp(name);
    load([pname]);
    i = pick;
    ei = para.ei(i);
    dir = [name,'-',theme,'-',model];
    name = dir;
%     loadData = false;
    testEE = true;
    testII = true;
    testEI = true;
    testIE = true;
%     multipleInput = false;
    multipleInput = true;
%     simpleTest = true;
    simpleTest = false;
    test = false;
%    test = true;
    fE
    fI
    seed = 132423;
    ignore = 0;
    ignorefE = 0;
    ignorefI = 0; 
    rateE = rateE/1000; % Hz/1000
    rateI = rateI/1000;
    para.fCurrent = 0;
    v0id = find(abs(v0 - 0.0)<1e-14);
    if isempty(v0id)
        disp('need vrest in vrange');
        assert(~isempty(v0id));
    end
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
            poolobj = gcp('nocreate');
            disp(poolobj.NumWorkers);
            feature('numcores');
        end
    else
        delete(gcp('nocreate'));
    end
    assert(nv0>1);
    nE = length(fE);
    nI = length(fI);
    % 1 2 3
    % E I EI
    diri = [dir,'/',num2str(i)];
    vRange = para.vRest(i) + (para.vT(i)-para.vRest(i))*v0;
    assert(vRange(v0id) == para.vRest(i));
    if ~loadData
        if exist(['./',dir])~=7
            disp(['making directory ',dir]);
            mkdir(dir);
        end
        if ~exist(['./',diri],'dir')
            disp(['making directory ',diri]);
            mkdir(diri);
        end
        if ~singleStored
            tic;
            if loadExplore
                load(['../library/explore-',name,'-',num2str(i),'th.mat'],'dtRange','dur');
            end
            [sEPSP,sIPSP,E_tmax,I_tmax,vleakage,dtRange,dur,sVcapE,sVcapI] = sPSP_check(silico,tstep,vRange,fI',fE',para,bool,name,dur,i,v0id,dtRange);
            dur0 = dur; % linear length
            t = 0:tstep:dur;
            t0 = 0:tstep:dur0;
            nt = round(dur/tstep)+1;
            nt0 = round(dur0/tstep)+1;
            ndt = length(dtRange);
            idtRange = round(dtRange/tstep)+1;
            assert(dtRange(ndt)<dur);
            for idt=2:ndt
                assert((dur - dtRange(idt))>=(dtRange(ndt-1)-dtRange(idt-1)));
            end
            toc;
            tic;
            sEPSP0 = squeeze(sEPSP(:,:,1,:));
            sIPSP0 = squeeze(sIPSP(:,:,1,:));
            vleakage0 = squeeze(vleakage(:,1,:));
            h = plotsPSP0(sEPSP0,sIPSP0,vleakage0,fE,fI,nv0,dur,nt,v0id,idtRange);
            if draw
                fname = ['sPSP0','-',name,'-',num2str(i)];
                printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
            end
            %% get nSPSP
            h = figure;
            maxsE = repmat(max(abs(squeeze(sEPSP(:,:,1,:)))),[nt,1]);
            maxsI = repmat(max(abs(squeeze(sIPSP(:,:,1,:)))),[nt,1]);
            for idt=1:ndt
                tpick = idtRange(idt):nt;
                sE = abs(squeeze(sEPSP(tpick,:,idt,:)));
                sI = abs(squeeze(sIPSP(tpick,:,idt,:)));
                nsEPSP = sE./maxsE(tpick,:,:);
                nsIPSP = sI./maxsI(tpick,:,:);
                subplot(2,2,1)
                hold on
                c1 = linspace(0,5/6,nv0);
                c2 = linspace(0.3,1,nE)';
                c3 = 0.5+0.3*(idt/ndt);
                for iF = 1:nE
                    for iv0 = 1:nv0
                        plot(t(tpick),squeeze(nsEPSP(:,iF,iv0)),'Color',hsv2rgb([c1(iv0),c2(iF),c3]));
                    end
                end
                title('norm. EPSP, hue v0, satur f');
                xlim([0,dur]);
                subplot(2,2,3)
                hold on
                for iF = 1:nE
                    for iv0 = 1:nv0
                        plot(dtRange(idt),log10(max(squeeze(sE(:,iF,iv0)))),'*','Color',hsv2rgb([c1(iv0),c2(iF),c3]));
                    end
                end
                ylabel('log10');
                subplot(2,2,2)
                hold on
                c2 = linspace(0.3,1,nI)';
                for iF = 1:nI
                    for iv0 = 1:nv0
                        plot(t(tpick),squeeze(nsIPSP(:,iF,iv0)),'Color',hsv2rgb([c1(iv0),c2(iF),c3]));
                    end
                end
                title('norm. IPSP, hue v0, satur f');
                xlim([0,dur]);
                subplot(2,2,4);
                hold on
                for iF = 1:nI
                    for iv0 = 1:nv0
                        plot(dtRange(idt),log10(max(squeeze(sI(:,iF,iv0)))),'*','Color',hsv2rgb([c1(iv0),c2(iF),c3]));
                    end
                end
                ylabel('log10');
            end
            if draw
                fname = ['nsPSP','-',name,'-',num2str(i)];
                printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
            end
            if explore
                save(['../library/explore-',name,'-',num2str(i),'th.mat'],'dtRange','dur');
                return;
            end
            toc;
            tic;
            [extraEPSP,extraIPSP,edtRange,nedt] = get_extraPSP(vleakage,silico,tstep,vRange,fI',fE',para,bool,name,dur,i,v0id,dtRange);
            toc;
            save(['../library/single-',name,'-',num2str(i),'th.mat'],'sEPSP','sIPSP','extraEPSP','extraIPSP','edtRange','nedt','vleakage','sEPSP0','sIPSP0','E_tmax','I_tmax','vleakage0','dtRange','idtRange','ndt','dur','dur0','nt','nt0','t','t0','-v7.3');
        else
            load(['../library/single-',name,'-',num2str(i),'th.mat']);
            h = plotsPSP0(sEPSP0,sIPSP0,vleakage0,fE,fI,nv0,dur,nt,v0id,idtRange);
        end
        dtRange
        edtRange
        disp(datestr(datetime('now')));
        idtCase = [1,round(ndt/2),ndt];
        ndtplot = length(idtCase);
        
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
        copyfile('noAdapV_k4.m',[dir,'/noAdapV-',name,'-',num2str(i),'th.m']);
        %   draw normalized EPSP
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
        if draw
            fname = ['nsPSP0','-',name,'-',num2str(i)];
            printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
        end
        
        % taking k at max |response|
%         [~,tQE] = max(sEPSP);
%         [~,tQI] = max(sIPSP);
%         tQE = round(mean(squeeze(tQE)));
%         tQI = round(mean(squeeze(tQI)));
        nF = nE + nI;
        kV  = zeros(nt,ndt,nF,nF,ndt,nv0);
        kEE  = zeros(nt,ndt,ndt,nv0);
        kII  = zeros(nt,ndt,ndt,nv0);
        kEI  = zeros(nt,ndt,ndt,nv0);
        kIE  = zeros(nt,ndt,ndt,nv0);

        if draw
            if testEE
                hEE0 = gobjects(1,ndtplot);
            end
            if testII
                hII0 = gobjects(1,ndtplot);
            end
            if testEI
                hEI0 = gobjects(1,ndtplot);
            end
            if testIE
                hIE0 = gobjects(1,ndtplot);
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
                [kV(:,:,1:nE,1:nE,idt,:),kEE(:,:,idt,:),pv,vadd,vEEDoublet,vEE,h,hh] = doubleCheck(silico,name,para,vv0,bool,tstep,iidt,dur,i,vleakage,sEPSP,sEPSP,extraEPSP,vRange,pp0,idt,ndt,dtRange,edtRange);
                if dtplot
                    hEE0(idtplot) = drawExample(pv,vadd,vEEDoublet,squeeze(vleakage(:,idt,:)),t,iidt,E_tmax(1,1,v0id),kEE(:,idt,idt,:),nE*nE,vRange,para.vRest(i),vEE,'V_{E}V_{E}','EE');
                end
                if pp0
                    if length(h) > 0
                        fname = [h.FileName,'-EE'];
                        printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
                    end
                    if length(hh) > 0
                        fname = [hh.FileName,'-EE'];
                        printpic(hh,diri,fname,picformat,printDriver,dpi,pos0);
                    end
                    fname = [num2str(idt),'k-EE'];
                    printpic(hEE0(idtplot),diri,fname,picformat,printDriver,dpi,pos0);
                end
            end
            if testII
                % II
                vv0 = vRange(v0id)*ones(nI*nI,1);
                para.tI = repmat([0,dtRange(idt)],[nI*nI,1]);
                para.f_I = [reshape(ones(nI,1)*fI,[nI*nI,1]),repmat(fI',[nI,1])];
                para.tE = zeros(size(para.tI));
                para.f_E = zeros(size(para.f_I));
                [kV(:,:,nE+1:nF,nE+1:nF,idt,:),kII(:,:,idt,:),pv,vadd,vIIDoublet,vII,h,hh] = doubleCheck(silico,name,para,vv0,bool,tstep,iidt,dur,i,vleakage,sIPSP,sIPSP,extraIPSP,vRange,pp0,idt,ndt,dtRange,edtRange);
                if dtplot
                    hII0(idtplot) = drawExample(pv,vadd,vIIDoublet,squeeze(vleakage(:,idt,:)),t,iidt,I_tmax(1,1,v0id),kII(:,idt,idt,:),nI*nI,vRange,para.vRest(i),vII,'V_{I}V_{I}','II');
                end
                if pp0
                    if length(h) > 0
                        fname = [h.FileName,'-II'];
                        printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
                    end
                    if length(hh) > 0
                        fname = [hh.FileName,'-II'];
                        printpic(hh,diri,fname,picformat,printDriver,dpi,pos0);
                    end
                    fname = [num2str(idt),'k-II'];
                    printpic(hII0(idtplot),diri,fname,picformat,printDriver,dpi,pos0);
                end
            end
            if testEI
                % EI 
                vv0 = vRange(v0id)*ones(nE*nI,1);
                para.f_E = reshape(ones(nI,1)*fE,[nE*nI,1]);
                para.tE = zeros(nE*nI,1);
                para.f_I = repmat(fI',[nE,1]);
                para.tI = dtRange(idt)*ones(nE*nI,1);
                [kV(:,:,nE+1:nF,1:nE,idt,:),kEI(:,:,idt,:),pv,vadd,vEIDoublet,vEI,h,hh] = doubleCheck(silico,name,para,vv0,bool,tstep,iidt,dur,i,vleakage,sEPSP,sIPSP,extraIPSP,vRange,pp0,idt,ndt,dtRange,edtRange);
                if dtplot
                    hEI0(idtplot) = drawExample(pv,vadd,vEIDoublet,squeeze(vleakage(:,idt,:)),t,iidt,I_tmax(1,1,v0id),kEI(:,idt,idt,:),nE*nI,vRange,para.vRest(i),vEI,'V_{E}V_{I}','EI');
                end
                if pp0
                    if length(h) > 0
                        fname = [h.FileName,'-EI'];
                        printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
                    end
                    if length(hh) > 0
                        fname = [hh.FileName,'-EI'];
                        printpic(hh,diri,fname,picformat,printDriver,dpi,pos0);
                    end
                    fname = [num2str(idt),'k-EI'];
                    printpic(hEI0(idtplot),diri,fname,picformat,printDriver,dpi,pos0);
                end
            end
            if testIE
                % IE
                vv0 = vRange(v0id)*ones(nE*nI,1);
                para.f_I = reshape(ones(nE,1)*fI,[nE*nI,1]);
                para.tI = zeros(nE*nI,1);
                para.f_E = repmat(fE',[nI,1]);
                para.tE = dtRange(idt)*ones(nE*nI,1);
                [kV(:,:,1:nE,nE+1:nF,idt,:),kIE(:,:,idt,:),pv,vadd,vIEDoublet,vIE,h,hh] = doubleCheck(silico,name,para,vv0,bool,tstep,iidt,dur,i,vleakage,sIPSP,sEPSP,extraEPSP,vRange,pp0,idt,ndt,dtRange,edtRange);
                if dtplot
                    hIE0(idtplot) = drawExample(pv,vadd,vIEDoublet,squeeze(vleakage(:,idt,:)),t,iidt,E_tmax(1,1,v0id),kIE(:,idt,idt,:),nE*nI,vRange,para.vRest(i),vIE,'V_{I}V_{E}','IE');
                end
                if pp0
                    if length(h) > 0
                        fname = [h.FileName,'-IE'];
                        printpic(h,diri,fname,picformat,printDriver,dpi,pos0);
                    end
                    if length(hh) > 0
                        fname = [hh.FileName,'-IE'];
                        printpic(hh,diri,fname,picformat,printDriver,dpi,pos0);
                    end
                    fname = [num2str(idt),'k-IE'];
                    printpic(hIE0(idtplot),diri,fname,picformat,printDriver,dpi,pos0);
                end
            end
        end
        toc;
        disp('kV generated');
        save(['../library/',name,'-',num2str(i),'th'],'kV','kEE','kIE','kEI','kII','vleakage','sEPSP','sIPSP','dur','vRange','fE','fI','nE','nI','ndt','i','dtRange','sEPSP0','sIPSP0','vleakage0','ei','tstep','dir','nt0','E_tmax','I_tmax','-v7.3');
        if draw
            assert(idtplot==ndtplot);
        end
    else
        if multipleInput
%             load(['single-',name,'-',num2str(i),'th']);
            load([name,'-',num2str(i),'th']);
            t0 = ((1:size(sEPSP,1))-1)*tstep;
            disp('kV loaded');
            dtRange
        end
    end
    if multipleInput
        disp('testing multiple inputs');
        % uniform K
            %>2 input
            % poisson inputs
        %  seed = 237587; % spike
%         poinumber = true;
        poinumber = true;
        poiend = false;
        idtRange = round(dtRange/tstep)+1;
        durpsp = dtRange(ndt-1) + dur - dtRange(ndt);
%         seed = 122435;
        rng(seed);
        l0 = round(durpsp/tstep);

        if ~simpleTest
%             dur = max(dur*1.5,2500); %ms
            dur = mdur;
            nXE = max(ceil(rateE*dur*2),10);
            nXI = max(ceil(rateI*dur*2),10);
            xE = rand(nXE,1);
            xI = rand(nXI,1);
            tE = zeros(size(xE));
            tI = zeros(size(xI));
            pfE = randi(nE-ignorefE,nXE,1);
            pfI = randi(nI-ignorefI,nXI,1);
            if rateE > 0
                t = 0;
                it = 1;
                t = t - log(xE(it))/rateE;
                while t < dur
                    tE(it) = t;
                    it = it + 1;
                    t = t - log(xE(it))/rateE;
                end
                it = it - 1;
                para.f_E = fE(pfE(1:it));
                nfE = it;
                para.tE = round(tE(1:it)'/tstep)*tstep;
            else
                para.tE = [];
                para.f_E = [];
            end
            if rateI > 0
                t = 0;
                it = 1;
                t = t - log(xI(it))/rateI;
                while t < dur
                    tI(it) = t;
                    it = it +1;
                    t = t - log(xI(it))/rateI;
                end
                it = it - 1;
                para.f_I = fI(pfI(1:it));
                nfI = it;
                para.tI = round(tI(1:it)'/tstep)*tstep;
            else
                para.tI = [];
                para.f_I = [];
            end
            fname = ['MultipleInputs',num2str(dur),'-E',num2str(rateE*1e3),'-I',num2str(rateI*1e3),'-S',num2str(seed),'-',name,'-',num2str(i),'th-v'];
        else
            pfE = [2];
            pfI = [2];
            nfE = length(pfE);
            nfI = length(pfI);
            para.f_E = fE(pfE);
            para.f_I = fI(pfI);
            para.tE = [10];
            para.tI = [5];
            dur = 600;
%             dur = 10;
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
        t = (0:length(v)-1)*tstep;
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
        if nt > size(vleakage0,1)
            vadd = [vleakage0(:,v0id);zeros(nt-size(vleakage0,1),1)+para.vRest(i)];
        else
            vadd = vleakage0(1:nt,v0id);
        end
%         vadd = vadd+para.vRest(i);
        vpred = vadd;
        iI = 1;
        finishI = false;
        for iE = 1:nfE
            while para.tI(iI) < para.tE(iE) && ~finishI
                if linear0
                    vI(:,iI) = sIPSP0(:,pfI(iI),v0id);
                else
                    vI(:,iI) = interpPSP(sIPSP0(:,pfI(iI),:),vadd(tIs(iI)),vRange);
                end
                tpick = tIs(iI):tIeadd(iI);
                vadd(tpick) = vadd(tpick) + vI(1:tIladd(iI),iI);
                if iI == nfI
                    finishI = true;
                    break;
                else
                    iI = iI + 1;
                end
            end
            if linear0
                vE(:,iE) = sEPSP0(:,pfE(iE),v0id);
            else
                vE(:,iE) = interpPSP(sEPSP0(:,pfE(iE),:),vadd(tEs(iE)),vRange);
            end
            tpick = tEs(iE):tEeadd(iE);
            vadd(tpick) = vadd(tpick) + vE(1:tEladd(iE),iE);
        end
        if ~finishI
            iI0 = iI;
            for iI = iI0:nfI
                if linear0
                    vI(:,iI) = sIPSP0(:,pfI(iI),v0id);
                else
                    vI(:,iI) = interpPSP(sIPSP0(:,pfI(iI),:),vadd(tIs(iI)),vRange);
                end
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

        finishI = false;
        iI = 1;
%         debug = true;
        debug = false;
        if debug
            dbstop if error
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
                        if round((para.tI(iI) - para.tE(iiE))/tstep) < idtRange(ndt) && round((para.tI(iI) - para.tE(iiE))/tstep)  >= 0
                            if debug && rule(iiE,0,0,iI)
                                hDebug = debugPreplot(hDebug,t,vpred);
                            end
                            [vpred, bEI{iiE,iI},tEI(:,iiE,iI)] = pred(vpred,iiE,iI,squeeze(sEPSP(:,pfE(iiE),:,:)),vI,...
                                                                      tEs,tEe,tEl,tIs,kEI,...
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
                        if round((para.tI(iI) - para.tI(iiI))/tstep) < idtRange(ndt)
                            if debug && rule(0,iiI,0,iI)
                                hDebug = debugPreplot(hDebug,t,vpred);
                            end
                            [vpred, bII{iiI,iI},tII(:,iiI,iI)] = pred(vpred,iiI,iI,squeeze(sIPSP(:,pfI(iiI),:,:)),vI,...
                                                                      tIs,tIe,tIl,tIs,kII,...
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
                for iiI = 1:iI
                    if round((para.tE(iE) - para.tI(iiI))/tstep) < idtRange(ndt) && round((para.tE(iE) - para.tI(iiI))/tstep) > 0
                        if debug && rule(0,iiI,iE,0)
                            hDebug = debugPreplot(hDebug,t,vpred);
                        end
                        [vpred, bIE{iiI,iE},tIE(:,iiI,iE)] = pred(vpred,iiI,iE,squeeze(sIPSP(:,pfI(iiI),:,:)),vE,...
                                                                  tIs,tIe,tIl,tEs,kIE,...
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
                    if round((para.tE(iE) - para.tE(iiE))/tstep) < idtRange(ndt)
                        if debug && rule(iiE,0,iE,0)
                            hDebug = debugPreplot(hDebug,t,vpred);
                        end
                        [vpred, bEE{iiE,iE},tEE(:,iiE,iE)] = pred(vpred,iiE,iE,squeeze(sEPSP(:,pfE(iiE),:,:)),vE,...
                                                                  tEs,tEe,tEl,tEs,kEE,...
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
                        if round((para.tI(iI) - para.tE(iiE))/tstep) < idtRange(ndt) && round((para.tI(iI) - para.tE(iiE))/tstep) >= 0
                            if debug && rule(iiE,0,0,iI)
                                hDebug = debugPreplot(hDebug,t,vpred);
                            end
                            [vpred,bEI{iiE,iI},tEI(:,iiE,iI)] = pred(vpred,iiE,iI,squeeze(sEPSP(:,pfE(iiE),:,:)),vI,...
                                                                     tEs,tEe,tEl,tIs,kEI,...
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
                        if round((para.tI(iI) - para.tI(iiI))/tstep) < idtRange(ndt)
                            if debug && rule(0,iiI,0,iI)
                                hDebug = debugPreplot(hDebug,t,vpred);
                            end
                            [vpred,bII{iiI,iI},tII(:,iiI,iI)] = pred(vpred,iiI,iI,squeeze(sIPSP(:,pfI(iiI),:,:)),vI,...
                                                                     tIs,tIe,tIl,tIs,kII,...
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
                if bool(para.type(i),4)
                    disp('iL in green');
                    iL = squeeze(tmp(6,:,:)).^2.*squeeze(tmp(7,:,:)).*para.gL(i).*(v-para.vCA);
                    plot(t,iL,abs(iL),'-g');
                end
                if bool(para.type(i),5)
                    disp('iT in yellow');
                    iT = squeeze(tmp(8,:,:)).^2.*squeeze(tmp(9,:,:)).*para.gT(i).*(v-para.vCA);
                    plot(t,iT,abs(iT),'-y');
                end

            end
            ylim(ax(2),[0,inf]);
            axes(ax(1));
            hold(ax(1),'on');
            ax(1).Color = 'none';
        end
        plot(t,vadd,'--m');
        plot(t,vpred,':m');
        maxv = max([max(v),max(vadd),max(vpred)]);
        minv = min([min(vadd),min(v),min(vpred)]);
        sl = maxv - minv;
        y = [minv-sl*0.2,maxv + sl*0.2];
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
        errorbar(1,sign(sum(vadd-v)) *mean(abs(vadd-v)),std(abs(vadd-v)),'*m');
        e=errorbar(2,sign(sum(vpred-v))*mean(abs(vpred-v)),std(abs(vpred-v)),'*b');
        set(ax3,'XTick',[1,2]);
        set(ax3,'XTickLabel',{'linear';'bilinear'});
        legend([e],{'k'});
        ylabel('error mV');
        xlim([0.5,2.5]);
        subplot(2,2,4);
        hold on
        plot(t,vadd-v,'--m');
        plot(t,vpred-v,':m');
        plot(t,zeros(size(t)),':k');
        xlim([t(1),t(end)]);
        ylim([max(min([min(vadd-v),min(vpred-v),0]),-3),min(max([max(vadd-v),max(vpred-v),0]),3)]);
        xlabel('ms');
        ylabel('err mV');
        
        printpic(hM,diri,fname,picformat,printDriver,dpi,pos0,false);
    end
end
function [EPSP,IPSP,dtRange,ndt] = get_extraPSP(vleakage,silico,tstep,vRange,fIRange,fERange,para,bool,name,dur,i,v0id,dtRange0)
    ndt0 = length(dtRange0);
    dtRange = [];
    for idt=1:ndt0
        for jdt=1:idt-1
            new_dt = dtRange0(idt) - dtRange0(jdt);
            numDupNew = new_dt - dtRange;
            numDupOld = new_dt - dtRange0;
            if sum(abs(numDupNew)<1e-14) == 0 && sum(abs(numDupOld)<1e-14) == 0
                dtRange = [dtRange, new_dt];
            end
        end
    end
    dtRange = sort(dtRange)
    ndt = length(dtRange);
    nv0 = length(vRange);
    nE = length(fERange);
    nI = length(fIRange);
    nt = round(dur/tstep)+1;
    EPSP = zeros(nt,nE,ndt,nv0);
    IPSP = zeros(nt,nI,ndt,nv0);
    idtRange = round(dtRange/tstep) + 1;
    if ndt > 0
        para.f_E = 0;
        para.tE = para.f_E;
        para.f_I = 0;
        para.tI = para.f_I;
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
                tmpLeak = [zeros(idtRange(idt)-1,1);vleakage(1:(nt-idtRange(idt)+1),1,iv0)];
                tmpEv = squeeze(tmpE(1,:,:)) - repmat(tmpLeak,[1,nE]);
                [~,E_tmax(:,idt,iv0)] = max(tmpEv);
                EPSP(:,:,idt,iv0) = tmpEv;
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
                tmpLeak = [zeros(idtRange(idt)-1,1);vleakage(1:(nt-idtRange(idt)+1),1,iv0)];
                tmpIv = squeeze(tmpI(1,:,:)) - repmat(tmpLeak,[1,nI]);
                [~,I_tmax(:,idt,iv0)] = min(tmpIv);
                IPSP(:,:,idt,iv0) = tmpIv;
            end
        end
    end
    disp('extra sPSP complete');
end
function [EPSP,IPSP,E_tmax,I_tmax,vleakage,dtRange,dur,sVcapE,sVcapI] = sPSP_check(silico,tstep,vRange,fIRange,fERange,para,bool,name,dur,i,v0id,dtRange0)
    nv0 = length(vRange);
    nE = length(fERange);
    nI = length(fIRange);
    thres0 = tstep^4;
    thres = thres0;
    if dur < 0
        difference = 1;
        dur = max(-dur,150);
        while true
            nt = round(dur/tstep)+1
            v0 = vRange(v0id);
            para.vtime = -1;
            para.f_E = max(fERange);
            para.tE = 0;
            para.tI = 0;
            para.f_I = 0; 
            tmpE = silico(name,v0,para,bool,tstep,dur,i,false);
            tmpE = abs(squeeze(tmpE(1,:,:))-v0);
            eratio = tmpE./max(tmpE);
            eratio(nt)
            if eratio(nt) < thres0;
                iE = nt;
                for ie=nt-1:-1:floor(nt/2)
                    if eratio(ie) > thres0
                        iE = ie+1;
                        break
                    end
                end
            else
                iE = [];
            end

            para.f_I = max(fIRange);
            para.tE = 0;
            para.tI = 0;
            para.f_E = 0; 
            tmpI = silico(name,v0,para,bool,tstep,dur,i,false);
            tmpI = abs(squeeze(tmpI(1,:,:))-v0);
            iratio = tmpE./max(tmpI);
            iratio(nt)
            if iratio(nt) < thres0;
                iI = nt;
                for ii=nt-1:-1:floor(nt/2)
                    if iratio(ii) > thres0
                        iI = ii+1;
                        break
                    end
                end
            else
                iI = [];
            end

            para.f_E = 0;
            para.tE = para.f_E;
            para.f_I = 0;
            para.tI = para.f_I;
            vleakTmp = zeros(nt,nv0);
            parfor iv0 = 1:nv0
                param = para;
                v0 = vRange(iv0);
                leakage = silico(name,v0,param,bool,tstep,dur,i,false);
                vleakTmp(:,iv0) = squeeze(leakage(1,:,:));
            end
            vltmp = abs(vleakTmp-vRange(v0id));
            lratio = vltmp./(ones(nt,1)*max(vltmp));
            iL = -1;
            for iv0 = 1:nv0
                if iv0 == v0id, continue; end
                lratio(nt,iv0)
                if lratio(nt,iv0) < thres0;
                    iLtmp = nt;
                    for il=nt-1:-1:floor(nt/2)
                        if lratio(il,iv0) > thres0
                            iLtmp = il+1;
                            break
                        end
                    end
                else
                    iLtmp = [];
                end
                if ~isempty(iLtmp)
                    if iLtmp > iL
                        iL = iLtmp;
                    end
                end
            end
            if ~isempty(iE) && ~isempty(iI) && iL~=-1
                nt = max([iE,iI,iL])
                vleakTmp = vleakTmp(1:nt,:);
                dur = (nt-1)*tstep;
                break;
            else
                dur = dur * 2;
            end
        end 
    else
        nt = round(dur/tstep)+1;
        para.f_E = 0;
        para.tE = para.f_E;
        para.f_I = 0;
        para.tI = para.f_I;
        para.vtime = -1;
        vleakTmp = zeros(nt,nv0);
        parfor iv0 = 1:nv0
            param = para;
            v0 = vRange(iv0);
            leakage = silico(name,v0,param,bool,tstep,dur,i,false);
            vleakTmp(:,iv0) = squeeze(leakage(1,:,:));
        end
    end
    disp('leakage complete');
    if dtRange0(1) < 0
        disp(['auto searching dtRange of length ', num2str(-dtRange0(1))]);
        t=0:tstep:dur;
        gE = G_E_df(0,t,fERange(nE),para.tau_er,para.tau_e,para.gConE);
        gI = G_I_df(0,t,fIRange(nI),para.tau_ir,para.tau_i,para.gConI);
        g = (gE'./max(gE)+gI'./max(gI))/2;

        ndt = -dtRange0(1);

        nnt = nt;
        nnt_old = 0;
        dp = 0.85;
        [idtRange, ~] = autoSpacing(g(1:nnt),ndt,round(nt/7),dp);
        qc = 0;
        while (true)
            ndt = length(idtRange)
            dtRange = idtRange*tstep
            vleakage = zeros(nt,ndt,nv0);
            E_tmax =  zeros(nE,ndt,nv0);
            I_tmax = zeros(nI,ndt,nv0);
            ampMaxE = zeros(nE,ndt,nv0);
            ampMaxI = zeros(nI,ndt,nv0);
            sVcapE = nv0 + ones(nE,ndt);
            sVcapI = nv0 + ones(nI,ndt);

            for idt = 1:ndt
                vleakage(:,idt,:) =[zeros(idtRange(idt),1,nv0);reshape(vleakTmp(1:nt-idtRange(idt),:),[nt-idtRange(idt),1,nv0])];
                %idtRange(idt)
                %if idtRange(idt) > 0
                %    disp(squeeze(vleakage(idtRange(idt),idt,:))');
                %end
                %disp(squeeze(vleakage(idtRange(idt)+1,idt,:))');
            end
            idtE = ndt;
            idtI = ndt;
            EPSP = zeros(nt,nE,ndt,nv0);
            IPSP = zeros(nt,nI,ndt,nv0);
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
                    [ampMaxE(:,idt,iv0),E_tmax(:,idt,iv0)] = max(abs(tmpEv(idtRange(idt)+1:nt,:)));
                    EPSP(:,:,idt,iv0) = tmpEv;
                end
                for idt = 1:ndt
                    [mTmp, iE] = max(ampMaxE(:,idt,iv0));
                    if mTmp>30
                        disp(mTmp);
                        if iv0 < sVcapE(iE,idt)
                            sVcapE(iE,idt) = iv0;
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
                    [ampMaxI(:,idt,iv0),I_tmax(:,idt,iv0)] = max(abs(tmpIv(idtRange(idt)+1:nt,:)));
                    IPSP(:,:,idt,iv0) = tmpIv;
                end
                for idt = 1:ndt
                    [mTmp, iI] = max(ampMaxI(:,idt,iv0));
                    if mTmp>30
                        disp(mTmp);
                        if iv0 < sVcapI(iI,idt)
                            sVcapI(iI,idt) = iv0;
                        end
                    end
                end
                ratio = ampMaxE(1,:,iv0)./ampMaxE(1,1,iv0);
                min(ratio)
                idtEtmp = sum(ratio<thres);
                if idtEtmp < idtE
                    idtE = idtEtmp;
                end
                ratio = ampMaxI(1,:,iv0)./ampMaxI(1,1,iv0);
                min(ratio)
                idtItmp = sum(ratio<thres);
                if idtItmp < idtI
                    idtI = idtItmp;
                end
            end
            if min(idtI,idtE) == 1 || qc == 5
                break;
            else
                if min(idtI,idtE) == 0 
                    idtRange(ndt+1) = idtRange(ndt) + (idtRange(ndt)-idtRange(ndt-1)) * 2;
                    qc = qc + 1
                    if nt - idtRange(ndt+1) < idtRange(ndt+1) - idtRange(ndt)
                        idtRange(ndt+1) = floor((nt + idtRange(ndt))/2)
                        if idtRange(ndt+1) == idtRange(ndt)
                            idtRange = idtRange(1:ndt);
                        end
                        qc = 5;
                        disp('cant expand any more');
                    end
                else
                    nnt = round((idtRange(min(idtI,idtE))+1 + nnt)/2)
                    [idtRange, ~] = autoSpacing(g(1:nnt),-dtRange0(1),round(nt/7),dp);
                    qc = qc + 1
                end
            end
        end
    else
        dtRange = dtRange0;
        idtRange = round(dtRange/tstep)
        ndt = length(dtRange);
        EPSP = zeros(nt,nE,ndt,nv0);
        E_tmax =  zeros(nE,ndt,nv0);
        IPSP = zeros(nt,nI,ndt,nv0);
        I_tmax = zeros(nI,ndt,nv0);
        vleakage = zeros(nt,ndt,nv0);
        ampMaxE = zeros(nE,ndt,nv0);
        ampMaxI = zeros(nI,ndt,nv0);
        sVcapE = nv0 + ones(nE,ndt);
        sVcapI = nv0 + ones(nI,ndt);
        for idt = 1:ndt
            vleakage(:,idt,:) =[zeros(idtRange(idt),1,nv0);reshape(vleakTmp(1:nt-idtRange(idt),:),[nt-idtRange(idt),1,nv0])];
            %if idtRange(idt) > 0
            %    disp(squeeze(vleakage(idtRange(idt),idt,:))');
            %end
            %disp(squeeze(vleakage(idtRange(idt)+1,idt,:))');
        end
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
                [ampMaxE(:,idt,iv0),E_tmax(:,idt,iv0)] = max(abs(tmpEv(idtRange(idt)+1:nt,:)));
                EPSP(:,:,idt,iv0) = tmpEv;
            end
            for idt = 1:ndt
                [mTmp, iE] = max(ampMaxE(:,idt,iv0));
                if mTmp>30
                    disp(mTmp);
                    if iv0 < sVcapE(iE,idt)
                        sVcapE(iE,idt) = iv0;
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
                [ampMaxI(:,idt,iv0),I_tmax(:,idt,iv0)] = max(abs(tmpIv(idtRange(idt)+1:nt,:)));
                IPSP(:,:,idt,iv0) = tmpIv;
            end
            for idt = 1:ndt
                [mTmp, iI] = max(ampMaxI(:,idt,iv0));
                if mTmp>30
                    disp(mTmp);
                    if iv0 < sVcapI(iI,idt)
                        sVcapI(iI,idt) = iv0;
                    end
                end
            end
        end
    end
    disp('sPSP complete');
end
function h = plotsPSP0(EPSP,IPSP,vleakage,fERange,fIRange,nv0,dur,nt,v0id,idtRange)
    nE = length(fERange);
    nI = length(fIRange);
    h = figure;
    s = [0.1,1.0];
    v = 1.0;
    t = linspace(0,dur,nt);
    subplot(2,2,2);
    plotSinglePSP(0,s,v,fERange,t,EPSP,nE,nv0,v0id,idtRange);
    xlabel('f');
    ylabel('t');
    zlabel('EPSP (mV)');
    subplot(2,2,4);
    plotSinglePSP(2/3,s,v,fIRange,t,IPSP,nI,nv0,v0id,idtRange);
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
function plotSinglePSP(h,s,v,x,y,z,m,n,j0,idtRange)
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
                plot3(zeros(length(idtRange),1)+x(i),y(idtRange),z(idtRange,i,j),':','Color','k','Marker','*');
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
            1
            2
            v1 = sPSP(:,:,1);
            v2 = sPSP(:,:,2);
            rv = (vTarget-vRange(2))/(vRange(1)-vRange(2))
            v = v2 + rv*(v1-v2);
        else
            ind-1
            ind
            v1 = sPSP(:,:,ind-1);
            v2 = sPSP(:,:,ind);
            rv = (vTarget-vRange(ind-1))/(vRange(ind)-vRange(ind-1))
            v = v1 + rv*(v2-v1);
        end
    end
end
function tmp = interpPSPdt(vec,i,tar,ref,l,eov)
    nref = length(ref);
    assert(tar<ref(nref));
    assert(tar>0);
    if tar < ref(1)
        i = 2;
        j = 1;
    else
        while ref(i+1) < tar
            i = i + 1;  
        end
        j = i+1;
    end
    r = (tar-ref(i))/(ref(j)-ref(i));
    base = vec(ref(i)+(0:(l-1)),:,i);
    tmp = base + r*(vec(ref(j)+(0:(l-1)),:,j)-base);
end
function [kV,k,pV,vaddV,vDoubletV,v1v2,h0,h00] = doubleCheck(silico,name,para,v0,bool,tstep,iidt,dur,i,vleakage,sPSP1,sPSP2,extraPSP2,vRange,pp,idt,ndt,dtRange,edtRange)
    global nvplot ndtplot iv0Case idtCase
    h0 = [];
    h00 = [];
    n1 = size(sPSP1,2);
    n2 = size(sPSP2,2);
    nt = round(dur/tstep)+1;
    nv0 = length(vRange);
    i2 = n2;
    i1 = n1;

    idtRange = round(dtRange/tstep)+1;
    nedt = length(edtRange);
    kV = zeros(nt,ndt,n2,n1,nv0);
    k = zeros(nt,ndt,nv0);
    r2 = zeros(nt,ndt,nv0);
    vDoubletV = zeros(nt,n1*n2,nv0);
    vaddV = zeros(n1*n2,nt,nv0);
    v1V = zeros(n1*n2,nt,nv0);
    v2V = zeros(n1*n2,nt,nv0);
    pV = zeros(n1*n2,nt,nv0);
    para.vtime = (iidt-1)*tstep;
    if idt == ndt
        dtOI = iidt:nt;
    else
        dtOI = iidt:(iidt+(nt-idtRange(idt+1)));
    end
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
           v2V(range,dtOI,iv0) = [sPSP2(1:dtl,:,1,iv0)]';
       end
       vaddV(:,dtOI,iv0) = v1V(:,dtOI,iv0) + v2V(:,dtOI,iv0);
       parfor it = dtOI
           [k(it,idt,iv0),pV(:,it,iv0),r2(it,idt,iv0)] = p_fit110k0(v1V(:,it,iv0),v2V(:,it,iv0),vDoubletV(it,:,iv0)');
           kV(it,idt,:,:,iv0) = reshape(vDoubletV(it,:,iv0)'-vaddV(:,it,iv0),[n2,n1]);
       end
    end

    v1v2 = v1V.*v2V;

    if pp && sum((idt-idtCase==0))==1
        h0 = figure;
        t = 0:tstep:dur;
        iv0 = nv0;
        range = (i1-1)*n2+(1:n2);
        subplot(2,1,1);
        h0.FileName = [num2str(idt), 'th dt'];
        title([num2str(iv0),'th v, ',h0.FileName]);
        pick = dtOI(1:round(length(dtOI)/2));
        difference = vDoubletV(pick,range(i2),iv0) - vaddV(range(i2),pick,iv0)';
        [ax,h1,h2] = plotyy(t(pick),vDoubletV(pick,range(i2),iv0),t(pick),difference);
        h1.Color = 'g';
        h1.LineStyle ='--';
        h1.LineWidth = 1.5;
        h2.Color = 'g';
        h2.LineStyle =':';
        h2.LineWidth = 1.5;
        ax(2).YColor = 'g';
        hold(ax(1));
        plot(ax(1),t(pick),v1V(range(1),pick,iv0),':r','LineWidth',1.5);
        plot(ax(1),t(pick),vaddV(range(i2),pick,iv0),':k','LineWidth',1.5);
        plot(ax(1),t(pick),pV(range(i2),pick,iv0),':m','LineWidth',1.5);
        plot(ax(1),t(pick),zeros(length(pick),1),'-k','LineWidth',1.5);
        plot(ax(1),t(pick),v2V(range(i2),pick,iv0),':b','LineWidth',1.5);
        ax(1).YLim = [-inf,inf];
        hold(ax(2));
        plot(ax(2),t(pick),zeros(length(pick),1),'-k');
        ax(1).XLim = [t(pick(1)),t(pick(end))];
        ax(2).XLim = [t(pick(1)),t(pick(end))];
        ylabel(ax(1),'PSP mV');
        ylabel(ax(2),'\Delta(PSP) mV');
        subplot(2,1,2);
        plot(t(pick),k(pick,idt,iv0));
        [ax,h1,h2] = plotyy(t(pick),k(pick,idt,iv0),t(pick),r2(pick,idt,iv0));
        ax(1).YColor = 'r';
        yl = zeros(1,2);
        yl(1) = min(k(pick,idt,iv0));
        yl(1) = max(-2.0,yl(1)-0.1*abs(yl(1)));
        yl(2) = max(k(pick,idt,iv0));
        yl(2) = min(2.0,yl(2)+0.1*abs(yl(2)));
        if yl(2) > yl(1)
            ax(1).YLim = yl;
        end
        ax(1).YTickMode = 'auto';
        ax(1).YLabel.String = 'k';
        ax(2).YColor = 'b';
        ax(2).YLim = [0,1.1];
        ax(2).YTickMode = 'auto';
        ax(2).YLabel.String = 'rsquare';
        h1.Color = 'r';
        h2.Color = 'b';
        ax(1).XLim = [t(pick(1)),t(pick(end))];
        ax(2).XLim = [t(pick(1)),t(pick(end))];
        xlabel('t ms');
    end
    vadd = zeros(n1*n2,nt,nv0);
    iedtRange = round(edtRange/tstep)+1;
    for jdt = (idt+1):ndt
        jjdt = idtRange(jdt);
        para.vtime = (jjdt-1)*tstep;
        if jdt == ndt
            iend = nt;
        else
            iend = jjdt + (nt - idtRange(jdt+1));
        end
        dtOI = jjdt:iend;
        dtl = length(dtOI);
        para.vtime = (jjdt-1)*tstep;
        vdoub = zeros(nt,n1*n2,nv0);
        parfor iv0 = 1:nv0
            param = para;
            param.newv = vRange(iv0);
            tmp = silico(name,v0,param,bool,tstep,dur,i,false);
            vdoub(:,:,iv0) = squeeze(tmp(1,:,:));
        end
        dt = dtRange(jdt) - dtRange(idt);
        pv = zeros(n1*n2,nt,nv0);
        for iv0 = 1:nv0
            vleak = vleakage(dtOI,jdt,iv0);
            vdoub(dtOI,:,iv0) = vdoub(dtOI,:,iv0) - repmat(vleak,[1,n1*n2]);
            if sum(abs(dt-dtRange)<1e-14)==1
                kdt = find(abs(dt-dtRange)<1e-14);
                assert(length(kdt)==1);
                kkdt = idtRange(kdt);
                tmp = squeeze(sPSP2(kkdt:(kkdt+dtl-1),:,kdt,iv0));
            else
                if sum(abs(dt-edtRange)<1e-14)==1
                    kdt = find(abs(dt-edtRange)<1e-14);
                    assert(length(kdt)==1);
                    kkdt = iedtRange(kdt);
                    tmp = squeeze(extraPSP2(kkdt:(kkdt+dtl-1),:,kdt,iv0));
                else
                    dt
                    edtRange
                    disp('should not reach here');
                    assert(false);
                end
            end
            for iF = 1:n1
                range = (iF-1)*n2+(1:n2);
                v1V(range,dtOI,iv0) = (sPSP1(dtOI,iF,jdt,iv0)*ones(1,n2))';
                v2V(range,dtOI,iv0) = tmp';
            end
            vadd(:,dtOI,iv0) = v1V(:,dtOI,iv0) + v2V(:,dtOI,iv0);
            parfor it = dtOI
                [k(it,jdt,iv0),pv(:,it,iv0),r2(it,jdt,iv0)] = p_fit110k0(v1V(:,it,iv0),v2V(:,it,iv0),vdoub(it,:,iv0)');
                kV(it,jdt,:,:,iv0) = reshape(vdoub(it,:,iv0)'-vadd(:,it,iv0),[n2,n1]);
            end
        end
        if pp && sum((jdt-idtCase==0))==1
            %disp('drawing dt here');
            h00 = figure;
            t = 0:tstep:dur;
            iv0 = nv0;
            range = (i1-1)*n2+(1:n2);
            subplot(2,1,1);
            h00.FileName = [num2str(idt),'-',num2str(jdt),'dt'];
            title([num2str(iv0),'th v, ',h00.FileName]);
            hold on
            pick = dtOI;
            plot(t(pick),v1V(range(1),pick,iv0),':r');
            plot(t(pick),v2V(range(i2),pick,iv0),':b');
            plot(t(pick),vadd(range(i2),pick,iv0),':k');
            plot(t(pick),pv(range(i2),pick,iv0),'--m');
            plot(t(pick),vdoub(pick,range(i2),iv0),'--g');
            plot(t(pick),zeros(length(pick),1),'-k');
            xlim([t(pick(1)),t(pick(end))]);
            ylabel('PSP mV');
            subplot(2,1,2);
            plot(t(pick),k(pick,jdt,iv0));
            [ax,h1,h2] = plotyy(t(pick),k(pick,jdt,iv0),t(pick),r2(pick,jdt,iv0));
            ax(1).YColor = 'r';
            yl = zeros(1,2);
            yl(1) = min(k(pick,jdt,iv0));
            yl(1) = max(-1.5,yl(1)-0.1*abs(yl(1)));
            yl(2) = max(k(pick,jdt,iv0));
            yl(2) = min(1.5,yl(2)+0.1*abs(yl(2)));
            if yl(2) > yl(1)
                ax(1).YLim = yl;
            end
            ax(1).YTickMode = 'auto';
            ax(1).YLabel.String = 'k';
            ax(2).YColor = 'b';
            ax(2).YLim = [0,1.1];
            ax(2).YTickMode = 'auto';
            ax(2).YLabel.String = 'rsquare';
            h1.Color = 'r';
            h2.Color = 'b';
            ax(1).XLim = [t(pick(1)),t(pick(end))];
            ax(2).XLim = [t(pick(1)),t(pick(end))];
            xlabel('t ms');
            pp = false;
        end
    end
end
function h = drawExample(pv0,vadd0,vDoublet0,vleakage0,t,iidt,tp,k0,n,vRange,vRest,v120,xl,titl)
    global nvplot iv0Case;
    nv0 = length(vRange);
    inv0 = 0;
    h = figure;
    nt = length(t);
    tpick=iidt:nt;
    tp = min(iidt + tp,nt);
    for iv0 = iv0Case;
        vleakage = repmat(vleakage0(:,iv0),[1,n]);
        k = k0(:,iv0);
        v12 = v120(:,:,iv0);
        vDoublet = vDoublet0(:,:,iv0) + vleakage - vRest;
        pv = pv0(:,:,iv0) + vleakage' - vRest;
        vadd = vadd0(:,:,iv0) + vleakage' - vRest;
        subplot(2,2,1);
        hold on
        for iF=1:n
            p = plot(t,vDoublet(:,iF));
            plot(t(tpick),pv(iF,tpick),':','Color',p.Color,'LineWidth',1.5);
            plot(t(tpick),vadd(iF,tpick),'--','Color',p.Color);
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
        if (inv0==length(iv0Case))
            ylabel('V_{SC}');
            xlabel(xl);
        end
        title({['k = ',num2str(k(tp),'%.3g'),' at v' num2str(iv0)]});
        inv0 = inv0 + 1;
    end
end
function printpic(h,dir,fname,picformat,printDriver,dpi,pos,closePic)
    if nargin < 8
        closePic = true;
    end
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
            disp(['saving fig:', fname]);
            print(h,[dir,'/',fname,'.',picformat],printDriver,'-loose',dpi);
        end
        saveas(h,[dir,'/',fname,'.fig']);
        if closePic
            close(h);
        end
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
    if isempty(fdt)
        tTarget
        idtRange
        assert(~isempty(fdt));
    end

    t0 = (idtRange(fdt)) + (0:lt);
    if fdt < ndt
        qdt = fdt + 1;
        t1 = (idtRange(qdt)) + (0:lt);
        rt = (tTarget-idtRange(fdt))/(idtRange(qdt)-idtRange(fdt))
    else
        qdt = fdt;
        t1 = t0;
        rt = 0
        assert(fdt == ndt);
    end
    if isempty(fv) % out of lower bound exterpolation
        k0 = k(t0,fdt,fdt,1);
        rv = (vTarget-v(1))/(v(1)-v(2))
        kVtmp = k0 + rv * (k0-k(t0,fdt,fdt,2)) + rt * (k(t1,qdt,qdt,1)-k0);
        sPSP0 = sPSP(t0,fdt,1);
        vtmp = sPSP0 + rv * (sPSP0-sPSP(t0,fdt,2)) + rt * (sPSP(t1,qdt,1)-sPSP0);
    else
        k0 = k(t0,fdt,fdt,fv);
        sPSP0 = sPSP(t0,fdt,fv);
        if fv == nv0  % out of upper bound exterpolation
            rv = (vTarget-v(fv))/(v(fv)-v(fv-1))
            kVtmp = k0 +  rv * (k0-k(t0,fdt,fdt,fv-1));
            vtmp = sPSP0 + rv * (sPSP0-sPSP(t0,fdt,fv-1));
        else
            rv = (vTarget-v(fv))/(v(fv+1)-v(fv))
            kVtmp = k0 + rv * (k(t0,fdt,fdt,fv+1)-k0);
            vtmp = sPSP0 + rv * (sPSP(t0,fdt,fv+1)-sPSP0);
        end
        kVtmp = kVtmp + rt * (k(t1,qdt,qdt,fv)-k0);
        vtmp = vtmp + rt * (sPSP(t1,qdt,fv)-sPSP0);
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
