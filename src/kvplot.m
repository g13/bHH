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
    diri = '.';
    load(name);
    nv0 = length(vRange);
    ndt = length(dtRange);
    i = j;
    nnv0 = 3;
    c = {'r','g','b'};
    h = figure;
    inv0 = 1;
    for iv0 = [1,round(nv0/2),nv0];
        subplot(1,nnv0,inv0);
        plot(dtRange,squeeze(kVEI(:,:,iv0))','.','Color',c{inv0});
        inv0 = inv0 + 1;
    end
    fname = ['kv-dt-',name,'-',num2str(i)];
    printpic_ext(h,diri,fname,picformat,printDriver,dpi,pos0);
    h = figure;
    nndt = 3;
    indt = 1;
    c = {'r','g','b'};
    for idt = [1,round(ndt/2),ndt-1];
        subplot(1,nnv0,indt);
        plot(vRange,squeeze(kVEI(:,idt,:))','.','Color',c{indt});
        indt = indt + 1;
    end
    fname = ['kv-v0-',name,'-',num2str(i)];
    printpic_ext(h,diri,fname,picformat,printDriver,dpi,pos0);
end
