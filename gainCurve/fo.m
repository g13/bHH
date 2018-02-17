function fo(picformat,cfgFn)
    if nargin < 2
        cfgFn = 'test.cfg';
        if nargin < 1
            picformat = 'fig';
        end
    end
    dir = '../';
    %  paper
    width = 11;     height = width/16*9;
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
    p = read_cfg(cfgFn);
    parts = strsplit(pwd, '/');
    parentfdr = parts{end};
    outputName = [p.theme,'-',parentfdr];
    ith = 1;
    poi_number = true;
    poi_end = true;
    diff_poi = false;
    load(p.lib_file,'dtRange','vRange','dur','kEI','sEPSP','fE','fI','nE','ndt');
    nbt = size(kEI,1);
    nv = length(vRange);
    ndur = size(sEPSP,1);
    edur = dtRange(ndt) - p.ignore_t;
    tstep = dur/(ndur-1);
    run_nt = round(p.run_t/tstep) + 1;
    DataFn = ['Data-',p.theme,'-',num2str(p.seed),'.bin'];
    tInFn = ['tIn-',p.theme,'-',num2str(p.seed),'.bin'];
    RasterFn = ['Raster-',p.theme,'-',num2str(p.seed),'.bin'];
    cpuFn = ['cpuTime-',p.theme,'-',num2str(p.seed),'.bin'];
    DataFid = fopen(DataFn);
    tmp = fread(DataFid,[1,inf],'double');
    size(tmp)
    clear tmp;
    fclose(DataFid);
    DataFid = fopen(DataFn);
    tInFid = fopen(tInFn);
    RasterFid = fopen(RasterFn);
    cpuFid = fopen(cpuFn);
    FontSize = 16;
    set(0,'DefaultAxesFontSize',FontSize);
    set(0,'DefaultTextFontSize',FontSize-2);
    fStr = [fE,fI];
    t = 0:tstep:p.run_t;
    rEl = length(p.rE);
    legendTable = {'sim','bi','li'};
    cpuTime = zeros(rEl,3);
    tspSize = zeros(rEl,3);
    rasterData = cell(rEl,3);
    err = zeros(rEl,2,2);
    for j = 1:rEl
        outputMat = fread(DataFid,[run_nt,8],'double');
        tspSize(j,1) = fread(RasterFid,[1,1], 'int');
        rasterData{j,1} = fread(RasterFid,[tspSize(j,1),1], 'double');
        tspSize(j,2) = fread(RasterFid,[1,1], 'int');
        rasterData{j,2} = fread(RasterFid,[tspSize(j,2),1], 'double');
        tspSize(j,3) = fread(RasterFid,[1,1], 'int');
        rasterData{j,3} = fread(RasterFid,[tspSize(j,3),1], 'double');
        cpuTime(j,1) = fread(cpuFid,[1,1],'double');
        cpuTime(j,2) = fread(cpuFid,[1,1],'double');
        cpuTime(j,3) = fread(cpuFid,[1,1],'double');
        simV = outputMat(:,1);
        biV  = outputMat(:,2);
        liV  = outputMat(:,3);
        gE   = outputMat(:,4);
        gI   = outputMat(:,5);
        ntmp = fread(tInFid,[1,1],'int');
        tin = fread(tInFid,[ntmp,1],'double');
        tID = fread(tInFid,[ntmp,1],'int64');
        tID = tID + 1;
        pE = (tID <= nE);
        pI = (tID > nE);
        tE = tin(pE);
        tI = tin(pI);
        Eid = tID(pE);
        Iid = tID(pI)-nE;
        Ein = length(tE);
        Iin = length(tI);
        figure;
        textFontSize = 8;
        minV = max([-80,min([min(simV),min(biV),min(liV)])]);
        maxV = min([-40,max([max(simV),max(biV),max(liV)])]);
        subplot(2,1,1);
        hold on
        plot(t,[simV,biV,liV]);
        xlabel('t ms');
        ylabel('mem mV');
        hold on
        vEtar = zeros(Ein,1);
        plot(tE,minV*ones(1,Ein),'.r');
        for i=1:Ein
            iv = floor(tE(i)/tstep)+1;
            jv = iv + 1;
            vEtar = simV(iv) + mod(tE(i),tstep)/tstep * (simV(jv)-simV(iv));
            plot([tE(i),tE(i)],[minV,vEtar],':r');
            if poi_number
                text(tE(i),minV+(vEtar-minV)*0.3,num2str(Eid(i)),'Color','r','FontSize',textFontSize);
                text(tE(i),minV+(vEtar-minV)*0.1,num2str(i),'Color','r','FontSize',textFontSize);
            end
            if poi_end 
                iv = floor((tE(i)+edur)/tstep);
                if iv+1 < run_nt
                    jv = iv + 1;
                    vEtar = simV(iv) + mod(tE(i),tstep)/tstep * (simV(jv)-simV(iv));
                    plot([tE(i),tE(i)]+edur,[vEtar,maxV],':r');
                    if poi_number
                        text(tE(i)+edur,maxV-(maxV-vEtar)*0.3,num2str(Eid(i)),'Color','r','FontSize',textFontSize);
                        text(tE(i)+edur,maxV-(maxV-vEtar)*0.1,num2str(i),'Color','r','FontSize',textFontSize);
                    end
                end
            end
        end
        vItar = zeros(Iin,1);
        plot(tI,minV*ones(1,Iin),'.b');
        for i=1:Iin
            iv = floor(tI(i)/tstep)+1;
            jv = iv + 1;
            vItar(i) = simV(iv) + mod(tI(i),tstep)/tstep * (simV(jv)-simV(iv));
            plot([tI(i),tI(i)],[minV,vItar(i)],':b');
            if poi_number
                text(tI(i),minV+(vItar(i)-minV)*0.3,num2str(Iid(i)),'Color','b','FontSize',textFontSize);
                text(tI(i),minV+(vItar(i)-minV)*0.1,num2str(i),'Color','b','FontSize',textFontSize);
            end
                
            if poi_end
                iv = floor((tI(i)+edur)/tstep);
                if iv+1 < run_nt
                    jv = iv + 1;
                    vtar = simV(iv) + mod(tI(i),tstep)/tstep * (simV(jv)-simV(iv));
                    plot([tI(i),tI(i)]+edur,[vtar,maxV],':b');
                    if poi_number
                        text(tI(i)+edur,maxV-(maxV-vtar)*0.3,num2str(Iid(i)),'Color','b','FontSize',textFontSize);
                        text(tI(i)+edur,maxV-(maxV-vtar)*0.1,num2str(i),'Color','b','FontSize',textFontSize);
                    end
                end
            end
        end
        legend(legendTable,'Location','northwest');
        ylim([minV,maxV]);
        subplot(2,2,3);
        hold on;
        signLiV = sign(sum(liV-simV));
        err(j,1,1) = mean(abs(liV-simV))*signLiV;
        err(j,2,1) = std(abs(liV-simV));
        errorbar(1,err(j,1,1),err(j,2,1),'b');
        plot(1,err(j,1,1),'*b')
        signBiV = sign(sum(biV-simV));
        err(j,1,2) = mean(abs(biV-simV))*signBiV;
        err(j,2,2) = std(abs(biV-simV));
        errorbar(2,err(j,1,2),err(j,2,2),'r');
        plot(2,err(j,1,2),'*r');
        xlim([0,3]);
        set(gca,'xtick',[1,2]);
        plot([0,3],[0,0],':k');
        set(gca,'xtick',[1,2]);
        set(gca,'xticklabel',{'linear','bilinear'});
        ylabel('error per timestep mV');
        subplot(2,2,4);
        dliV = liV - simV;
        sum(dliV)
        plot(t,dliV);
        hold on
        dbiV = biV - simV;
        sum(dbiV)
        plot(t,dbiV);
        plot(t,zeros(size(t)),':k');
        legend({'\Delta(li-sim)','\Delta(bi-sim)'});
        minV = min(min(dbiV),min(dliV));
        maxV = max(max(dbiV),max(dliV));
        vEtar = zeros(Ein,1);
        plot(tE,minV*ones(1,Ein),'.r');
        targetV = dliV;
        if diff_poi
            for i=1:Ein
                iv = floor(tE(i)/tstep)+1;
                jv = iv + 1;
                vEtar = targetV(iv) + mod(tE(i),tstep)/tstep * (targetV(jv)-targetV(iv));
                plot([tE(i),tE(i)],[minV,vEtar],':r');
                if poi_number
                    text(tE(i),minV+(vEtar-minV)*0.3,num2str(Eid(i)),'Color','r','FontSize',textFontSize);
                    text(tE(i),minV+(vEtar-minV)*0.1,num2str(i),'Color','r','FontSize',textFontSize);
                end
             
                if poi_end
                    iv = floor((tE(i)+edur)/tstep);
                    if iv+1 < run_nt
                        jv = iv + 1;
                        vEtar = targetV(iv) + mod(tE(i),tstep)/tstep * (targetV(jv)-targetV(iv));
                        plot([tE(i),tE(i)]+edur,[vEtar,maxV],':r');
                        if poi_number
                            text(tE(i)+edur,maxV-(maxV-vEtar)*0.3,num2str(Eid(i)),'Color','r','FontSize',textFontSize);
                            text(tE(i)+edur,maxV-(maxV-vEtar)*0.1,num2str(i),'Color','r','FontSize',textFontSize);
                        end
                    end
                end
            end
            vItar = zeros(Iin,1);
            plot(tI,minV*ones(1,Iin),'.b');
            for i=1:Iin
                iv = floor(tI(i)/tstep)+1;
                jv = iv + 1;
                vItar(i) = targetV(iv) + mod(tI(i),tstep)/tstep * (targetV(jv)-targetV(iv));
                plot([tI(i),tI(i)],[minV,vItar(i)],':b');
                if poi_number
                    text(tI(i),minV+(vItar(i)-minV)*0.3,num2str(Iid(i)),'Color','b','FontSize',textFontSize);
                    text(tI(i),minV+(vItar(i)-minV)*0.1,num2str(i),'Color','b','FontSize',textFontSize);
                end
                    
                if poi_end
                    iv = floor((tI(i)+edur)/tstep);
                    if iv+1 < run_nt
                        jv = iv + 1;
                        vtar = targetV(iv) + mod(tI(i),tstep)/tstep * (targetV(jv)-targetV(iv));
                        plot([tI(i),tI(i)]+edur,[vtar,maxV],':b');
                        if poi_number
                            text(tI(i)+edur,maxV-(maxV-vtar)*0.3,num2str(Iid(i)),'Color','b','FontSize',textFontSize);
                            text(tI(i)+edur,maxV-(maxV-vtar)*0.1,num2str(i),'Color','b','FontSize',textFontSize);
                        end
                    end
                end
            end
        end
        fname =[outputName,'-trial',num2str(j)];
        printpic(gcf,dir,fname,picformat,printDriver,dpi,pos0);
    end
    fclose(DataFid);
    fclose(tInFid);
    fclose(RasterFid);
    fclose(cpuFid);
    figure;
    subplot(1,2,1)
    plot(p.rE+p.rI,cpuTime);
    legend(legendTable);
    xlabel('E+I input rate');
    ylabel('cpu time cost');
    subplot(1,2,2)
    hold on
    errorbar(p.rE+p.rI,err(:,1,1),err(:,2,1),'b');
    errorbar((p.rE+p.rI)+0.2,err(:,1,2),err(:,2,2),'r');
    xlabel('E+I input rate');
    ylabel('error per timestep');
    legend({'linear','bilinear'});
    fname =[outputName,'-err_rate'];
    printpic(gcf,dir,fname,picformat,printDriver,dpi,pos0);
    save([outputName,'-Raster.mat'],'rasterData');
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
        close(h);
    end
end
