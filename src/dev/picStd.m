function [printDriver,dpi,pos0,picformat] = picStd(picformat,width,height,mleft,mbot,mright,mtop)
    if nargin < 2
        width = 17;
        if nargin < 3
            height = width/16*9;            
            if nargin < 4
                mleft = 0.0;
                if nargin < 5
                    mbot = 0.0;
                    if nargin < 6
                        mright = 0.0;
                        if nargin < 7
                            mtop = 0.0;
                        end
                    end
                end
            end
        end
    end

    if ~isempty(picformat)
        if strcmp(picformat,'psc2')
            printDriver = ['-de',picformat];
            picformat = 'eps';
        else
            printDriver = ['-d',picformat];
        end
        dpi = '-r300';
        %  paper
%         width = 17;     height = width/16*9;
        % figure
        left = mleft*width;
        bot = mbot*height;    
        fwidth = width-mleft-width*mright;
        fheight = height-mbot-height*mtop;
        pos0 = [width, height, left,bot,fwidth,fheight];
    else
        dpi = '';
        printDriver = '';
        pos0 = 0;
    end        
end