function kplot(name,picformat,j)
    if nargin < 2
        picformat = '';
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
    m = 0.5;
    load(name);
    i=j;
    diri = '.';
    ndt = length(dtRange);
    for idt = [1,round(ndt/2),ndt-1];
        h = kvplot(kVEI,dtRange,vRange,idt,m);
        fname = ['kVEI-',name,'-',num2str(i),'-',num2str(idt)];
        printpic_ext(h,diri,fname,picformat,printDriver,dpi,pos0);
        h = kqvplot(kQVEI,dtRange,vRange,fE,idt,m);
        fname = ['kQVEI-',name,'-',num2str(i),'-',num2str(idt)];
        printpic_ext(h,diri,fname,picformat,printDriver,dpi,pos0);
    end
end
% kVEI  = zeros(nt,ndt,nv0,nv0);
% kQVEI = zeros(nt,ndt,nv0,nE,nv0);
function h = kvplot(kv,dtRange,vRange,idt,m)
    nt = size(kv,1);
    ndt = length(dtRange);
    nv0 = length(vRange);
    assert(ndt == size(kv,2));
    assert(nv0 == size(kv,3));
    h = figure;
    subplot(2,2,1)
    hold on
    iv0 = 1;
    for iiv0 = 1:nv0
        plot(1:nt,squeeze(kv(:,idt,iv0,iiv0)),'k','LineWidth',iiv0*m);
    end
    title(['diff 2v0, sample 1v0 = ',num2str(iv0)]);
    subplot(2,2,3)
    hold on
    for iiv0 = 1:nv0
        plot(1:nt,std(squeeze(kv(:,idt,:,iiv0)),1,2),'k','LineWidth',iiv0*m);
    end
    ylabel('std over 1v0');

    %subplot(2,2,2)
    %hold on
    %for iv0 = 1:nv0
    %    plot(1:nt,squeeze(kv(:,idt,:,iv0)),'k','LineWidth',iv0*m);
    %end
    %title('2nd v0');
end
function h = kqvplot(kqv,dtRange,vRange,f,idt,m);
    nt = size(kqv,1);
    ndt = length(dtRange);
    nv0 = length(vRange);
    nf = length(f);
    assert(ndt == size(kqv,2));
    assert(nv0 == size(kqv,3));
    assert(nf == size(kqv,4));
    h = figure;
    hold on
    c = {'r','m','g','b'};
    for iF = 1:nf
        subplot(2,2,1)
        hold on
        iv0 = 1;
        for iiv0 = 1:nv0
            plot(1:nt,squeeze(kqv(:,idt,iv0,iF,iiv0)),c{iF},'LineWidth',iiv0*m);
        end
        title(['diff 2v0, sample 1v0 = ',num2str(iv0)]);
        subplot(2,2,3)
        hold on
        for iiv0 = 1:nv0
            plot(1:nt,std(squeeze(kqv(:,idt,:,iF,iiv0)),1,2),c{iF},'LineWidth',iiv0*m);
        end
        ylabel('std over 1v0');

        title('1st v0');
    end
    subplot(2,2,2)
    hold on
    iv0 = 5;
    for iiv0 = 1:nv0
        target = squeeze(kqv(:,idt,iv0,:,iiv0));
        target = std(target,1,2)./mean(target,2);
        target(isnan(target)) = 0;
        plot(1:nt,target,'k','LineWidth',iiv0*m);
    end
    ylabel('std/mean f');
    title(['diff 2nd v0, sample 1v0 = ',num2str(iv0)]);
end
