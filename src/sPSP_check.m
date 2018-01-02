function sPSP_check(name,picformat)
    width = 17;     height = width/16*9;
    % figure
    mleft = 0.0;    left = mleft*width;
    mbot = 0.0;     bot = mbot*height;    
    mright = 0.0;   fwidth = width-mleft-width*mright;
    mtop = 0.0;     fheight = height-mbot-height*mtop;
    pos = [width, height, left,bot,fwidth,fheight];  
    addpath('channels/');    
    if ~isempty(picformat)
        if strcmp(picformat,'psc2')
            printDriver = ['-de',picformat];
            picformat = 'eps';
        else
            printDriver = ['-d',picformat];
        end
        dpi = '-r300';
    end
    
    load(['parameters-',name,'.mat']);
    i = 1;
    tstep = 0.05;
    fI = [0, (0.25:0.25:1.0)*1e-6];
    fE = [0, (0.25:0.25:1.0)*1e-6];
    nE = length(fE);
    l0 = 250;
    dur = 2*l0;
    t = 0:tstep:dur;
    nt = round(dur/tstep)+1;
    %plotIndivid = true;
    plotIndivid = false;
    nSeq = nE*(nE+1)/2-1;
    iSeq = zeros(nSeq,2);
    iiE = 1;
    for iE = 1:nE
        for jE = iE:nE
            if iiE~=1|| jE ~=iE 
                iSeq(iiE,:) = [iE,jE];
                iiE = iiE+1;
            end
        end
    end
    assert(iiE-1==nSeq);
    tSeq = rand(nSeq,1)*10+15;
    tSeq(1:nE-1) = 0;
    ratio = zeros(nt,5,nSeq);
    ratio0 = ratio;
    tLeak = cell(nSeq,1);
    hR = figure;
    para.fCurrent = 0;
    para.f_I = 0;
    para.tI = 0;
    satur = linspace(0.3,1,nSeq);
    gv0.exist = true;
    for iE = 1:nSeq;
        disp(iE);
        para.tE = [0,tSeq(iE)];
        para.f_E = fE(iSeq(iE,:));
        v0 = para.vRest(i);
        yE = RK4_df(name,v0,para,bool,tstep,dur,i,false);
        vE = yE(1,:,i);
        mE = yE(2,:,i);
        nE = yE(3,:,i);
        hE = yE(4,:,i);
        pE = yE(5,:,i);
        qE = yE(6,:,i);
        rE = yE(7,:,i);
        sE = yE(8,:,i);
        uE = yE(9,:,i);
        [vEmax,tEmax] = max(vE);
    
        para.f_E = 0;
        para.tE = 0;
        gv0.m0 = mE(tEmax);
        gv0.n0 = nE(tEmax);
        gv0.h0 = hE(tEmax);
        gv0.p0 = pE(tEmax);
        gv0.q0 = qE(tEmax);
        gv0.r0 = rE(tEmax);
        gv0.s0 = sE(tEmax);
        gv0.u0 = uE(tEmax);
        yLeak = RK4_df(name,vEmax,para,bool,tstep,dur-(tEmax-1)*tstep,i,false);
        vLeak = yLeak(1,:,i);
        mLeak = yLeak(2,:,i);
        nLeak = yLeak(3,:,i);
        hLeak = yLeak(4,:,i);
        pLeak = yLeak(5,:,i);
        qLeak = yLeak(6,:,i);
        rLeak = yLeak(7,:,i);
        sLeak = yLeak(8,:,i);
        uLeak = yLeak(9,:,i);
    
        yLeak0 = RK4_df(name,vEmax,para,bool,tstep,dur-(tEmax-1)*tstep,i,false,gv0);
        vLeak0 = yLeak0(1,:,i);
        mLeak0 = yLeak0(2,:,i);
        nLeak0 = yLeak0(3,:,i);
        hLeak0 = yLeak0(4,:,i);
        pLeak0 = yLeak0(5,:,i);
        qLeak0 = yLeak0(6,:,i);
        rLeak0 = yLeak0(7,:,i);
        sLeak0 = yLeak0(8,:,i);
        uLeak0 = yLeak0(9,:,i);
    
        tLeak{iE} = t(tEmax:nt);
        assert(mLeak0(1) == mE(tEmax));
        assert(nLeak0(1) == nE(tEmax));
        assert(hLeak0(1) == hE(tEmax));
        assert(pLeak0(1) == pE(tEmax));
        assert(qLeak0(1) == qE(tEmax));
        assert(rLeak0(1) == rE(tEmax));
        assert(sLeak0(1) == sE(tEmax));
        assert(uLeak0(1) == uE(tEmax) || isnan(uE(tEmax)) && isnan(uLeak0(1)));
        if plotIndivid
            h = figure;
            subplot(2,1,1)
            hold on
            plot(t,vE,'r');
            plot(tLeak{iE},vLeak,'-.k');
            plot(tLeak{iE},vLeak0,':k');
            plot(t,para.vRest(i)*ones(size(t)),':g');
            subplot(2,1,2);
            hold on
            pM = plot(t,mE);
            plot(tLeak{iE},mLeak0,':','Color',pM.Color);
            plot(tLeak{iE},mLeak,'-.','Color',pM.Color);
            pN = plot(t,nE);
            plot(tLeak{iE},nLeak0,':','Color',pN.Color);
            plot(tLeak{iE},nLeak,'-.','Color',pN.Color);
            pH = plot(t,hE);
            plot(tLeak{iE},hLeak0,':','Color',pH.Color);
            plot(tLeak{iE},hLeak,'-.','Color',pH.Color);
            pP = plot(t,pE);
            plot(tLeak{iE},pLeak0,':','Color',pP.Color);
            plot(tLeak{iE},pLeak,'-.','Color',pP.Color);
    
            legend([pM,pN,pH,pP],{'m','n','h','p'});
            plot(t,ones(size(t)),':g');
            ylim([0,max([mLeak,nLeak,hLeak,pLeak,mE,nE,hE,pE])*1.2]);
            xlabel('time ms');
            fname = [name,'-',num2str(i),'th-',num2str(iE)];
            printpic_ext(h,fname,picformat,printDriver,dpi,pos);
        end
        vMin = min([vE(tEmax:nt),vLeak]);
        vvE = vE-vMin;
        vvLeak = vLeak - vMin;
    
        target = vvLeak./vvE(tEmax:nt);
        ratio(tEmax:nt,1,iE) = target;
        target = mLeak./mE(tEmax:nt);
        ratio(tEmax:nt,2,iE) = target;
        target = nLeak./nE(tEmax:nt);
        ratio(tEmax:nt,3,iE) = target;
        target = hLeak./hE(tEmax:nt);
        ratio(tEmax:nt,4,iE) = target;
        target = pLeak./pE(tEmax:nt);
        ratio(tEmax:nt,5,iE) = target;
    
        vMin0 = min([vLeak0,vLeak]);
        vvLeak = vLeak-vMin0;
        vvLeak0 = vLeak0 - vMin0;
    
        target = vvLeak./vvLeak0;
        ratio0(tEmax:nt,1,iE) = target;
        target = mLeak./mLeak0;
        ratio0(tEmax:nt,2,iE) = target;
        target = nLeak./nLeak0;
        ratio0(tEmax:nt,3,iE) = target;
        target = hLeak./hLeak0;
        ratio0(tEmax:nt,4,iE) = target;
        target = pLeak./pLeak0;
        ratio0(tEmax:nt,5,iE) = target;
        figure(hR);
        subplot(2,2,1)
        hold on
        plot(tLeak{iE},ratio(tEmax:nt,1,iE),'k');
        subplot(2,2,3)
        hold on
        plot(tLeak{iE},ratio0(tEmax:nt,1,iE),'k');
        
        subplot(2,2,2)
        hold on
        if iE == 1
            pM = plot(tLeak{iE},ratio(tEmax:nt,2,iE));
            pN = plot(tLeak{iE},ratio(tEmax:nt,3,iE));
            pH = plot(tLeak{iE},ratio(tEmax:nt,4,iE));
            pP = plot(tLeak{iE},ratio(tEmax:nt,5,iE));
            mcolor = pM.Color;
            mhsv = rgb2hsv(mcolor);
            mhsv(2) = satur(iE);
            pM.Color = hsv2rgb(mhsv);
            ncolor = pN.Color;
            nhsv = rgb2hsv(ncolor);
            nhsv(2) = satur(iE);
            pN.Color = hsv2rgb(nhsv);
            hcolor = pH.Color;
            hhsv = rgb2hsv(hcolor);
            hhsv(2) = satur(iE);
            pH.Color = hsv2rgb(hhsv);
            pcolor = pP.Color;
            phsv = rgb2hsv(pcolor);
            phsv(2) = satur(iE);
            pP.Color = hsv2rgb(phsv);
            legend([pM,pN,pH,pP],{'m','n','h','p'});
        else
            mhsv(2) = satur(iE);
            nhsv(2) = satur(iE);
            hhsv(2) = satur(iE);
            phsv(2) = satur(iE);
            plot(tLeak{iE},ratio(tEmax:nt,2,iE),'Color',hsv2rgb(mhsv));
            plot(tLeak{iE},ratio(tEmax:nt,3,iE),'Color',hsv2rgb(nhsv));
            plot(tLeak{iE},ratio(tEmax:nt,4,iE),'Color',hsv2rgb(hhsv));
            plot(tLeak{iE},ratio(tEmax:nt,5,iE),'Color',hsv2rgb(phsv));
        end
        subplot(2,2,4)
        hold on
        plot(tLeak{iE},ratio0(tEmax:nt,2,iE),'-.','Color',hsv2rgb(mhsv));
        plot(tLeak{iE},ratio0(tEmax:nt,3,iE),'-.','Color',hsv2rgb(nhsv));
        plot(tLeak{iE},ratio0(tEmax:nt,4,iE),'-.','Color',hsv2rgb(hhsv));
        plot(tLeak{iE},ratio0(tEmax:nt,5,iE),'-.','Color',hsv2rgb(phsv));
    end
    figure(hR)
    subplot(2,2,2)
    plot(t,ones(size(t)),':k');
    xlabel('time');
    ylabel('ratio leak/input');
    
    subplot(2,2,4)
    xlabel('time');
    ylabel('ratio leak0/leak');
    
    fname = [name,'-',num2str(i),'th-ratio'];
    printpic_ext(hR,fname,picformat,printDriver,dpi,pos);
end
