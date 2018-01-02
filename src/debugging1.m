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
plot(kV
diri = 'RS_exc_Rat-noAdap-HH/1/';
fname = ['testK-',num2str(iv0),'-',num2str(idt),'-',num2str(iiv0)];
