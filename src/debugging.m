h00 = figure;
picformat = 'png';
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
ii = 1;
range = (ii-1)*n2+(1:n2);
i2 = n2;
t = 0:tstep:dur;
subplot(2,2,1)
    hold on
    plot(t,v1V(range(1),:,iiv0),':r');
    plot(t,v2V(range(i2),:,iiv0),':b');
    plot(t,vaddV(range(i2),:,iiv0),':k');
    plot(t,vDoubletV(:,range(i2),iiv0),'k');
subplot(2,2,2)
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
    plot(ax(2),t(pick),kQV(ii,pick,iiv0),'--g');
    plot(ax(2),t(iQV(ii,iiv0)),kQV(ii,iQV(ii,iiv0),iiv0),'*g');
    ylim(ax(2),[-100,100]);
subplot(2,2,3)
    hold on
    pick2 = 1:nt-iidt+1;
    plot(t(pick),vt1V(pick2,range(i2),iiv0),':r');
    plot(t(pick),vt2V(pick2,range(i2),iiv0),':b');
    plot(t(pick),interpPSP2V(pick2,range(i2),iiv0),':k');
subplot(2,2,4)
    hold on
    tp = min(iidt + 22/tstep,nt);
    v12 = v1V(:,tp,iiv0).*v2V(:,tp,iiv0);
    p = plot(v12,(vDoubletV(tp,:,iiv0)'-vaddV(:,tp,iiv0)),'o');
    xx = sort(v12);
    plot(xx,xx*kV(tp,iiv0),':k','LineWidth',1.5);
    plot(xx,(kQV(:,tp,iiv0)*xx')','--','Color',p.Color);
    
diri = 'RS_exc_Rat-noAdap-HH/1/';
fname = ['test-',num2str(iv0),'-',num2str(idt),'-',num2str(iiv0)];
printpic_ext(h00,diri,fname,picformat,printDriver,dpi,pos0);
