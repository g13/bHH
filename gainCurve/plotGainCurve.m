function plotGainCurve(picformat,cfgFn,ld)
    if nargin < 3
        ld = true;
        if nargin < 2
            cfgFn = 'test.cfg';
            if nargin < 1
                picformat = 'fig';
            end
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
    rEl = length(p.rE);
    load(p.lib_file,'dtRange','vRange','dur','kEI','sEPSP','fE','fI','nE','ndt');
    nbt = size(kEI,1);
    nv = length(vRange);
    ndur = size(sEPSP,1);
    edur = dtRange(ndt) - p.ignore_t;
    tstep = dur/(ndur-1);
    run_nt = round(p.run_t/tstep) + 1;
    parts = strsplit(pwd, '/');
    parentfdr = parts{end};
    outputName = [p.theme,'-',parentfdr];
    if ld
        if exist([outputName,'-Raster.mat'])  
            load([outputName,'-Raster.mat']);
        else
            disp('no raster data exist');
        end
    end
    if ~exist([dir,outputName,'-Raster.mat']) || ~ld
        disp('collecting raster data');
        RasterFn = ['Raster-',p.theme,'-',num2str(p.seed),'.bin'];
        RasterFid = fopen(RasterFn);
        rasterData = cell(rEl,3);
        for j=1:rEl
            tspSize(j,1) = fread(RasterFid,[1,1], 'int');
            rasterData{j,1} = fread(RasterFid,[tspSize(j,1),1], 'double');
            tspSize(j,2) = fread(RasterFid,[1,1], 'int');
            rasterData{j,2} = fread(RasterFid,[tspSize(j,2),1], 'double');
            tspSize(j,3) = fread(RasterFid,[1,1], 'int');
            rasterData{j,3} = fread(RasterFid,[tspSize(j,3),1], 'double');
        end
        fclose(RasterFid);
        save([dir,outputName,'-Raster.mat'],'rasterData','tspSize');
    end
    for j=1:rEl
        subplot(rEl,2,2*(j-1)+1);
        hold on
        plot(rasterData{j,1},zeros(tspSize(j,1),1)+1,'.','MarkerSize',10);
        plot(rasterData{j,2},zeros(tspSize(j,2),1)+2,'.','MarkerSize',10);
        plot(rasterData{j,3},zeros(tspSize(j,3),1)+3,'.','MarkerSize',10);
        set(gca,'YTick',[1,2,3],'YTickLabel',{'sim','bi','li'});
        xlim([0,p.run_t]);
        ylim([0,4]);
    end
    subplot(1,2,2);
    hold on
    plot(p.rE+p.rI,tspSize(:,1)/p.run_t,'*');
    plot(p.rE+p.rI,tspSize(:,2)/p.run_t,'o');
    plot(p.rE+p.rI,tspSize(:,3)/p.run_t,'s');
    legend({'sim','bi','li'});
    fname =[outputName,'-active'];
    printpic(gcf,dir,fname,picformat,printDriver,dpi,pos0);
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
