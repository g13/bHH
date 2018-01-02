function [y,spiked] = RK2_EIF_df(name,v0,para,bool,tstep,etime,ipick,savedata,reset)  %RK4 method
    if nargin < 9
        reset = true;
        if nargin < 8
            savedata = false;
        end
    end
    [y,spiked] = RK2_IF_df(name,v0,para,bool,tstep,etime,ipick,savedata,true,reset);
end
