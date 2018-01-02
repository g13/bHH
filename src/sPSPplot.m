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
iF = 4;
nt = size(sEPSP,1);
v = -0.4:0.1:1.0;
nv0 = length(v);
idt = 8;
tstep = 0.1;
dtstep = round(280/40);
t = 0:0.1:280;
iidt = (idt-1)*(dtstep/tstep)+1;
assert(nv0==size(sEPSP,3));
hold on
c = [0,0,1.0];
s = linspace(0.1,1.0,nv0);
for iv0 = 1:nv0
    for iiv0 = 1:nv0
        c(2) = s(iiv0);
        z = squeeze(sEPSP(iidt:nt,iF,iv0,idt,iiv0));
        plot3(t(iidt:nt),v(iv0)*ones(nt-iidt+1,1),z,'Color',hsv2rgb(c));
    end
end
[mz, ind] = max(abs(z(:)));
campos([0.1*max(t),0.2*max(v),0.8*sign(z(ind))*mz]);

fname = ['single-sEPSP-',num2str(idt)];
printpic_ext(h00,'.',fname,picformat,printDriver,dpi,pos0);
