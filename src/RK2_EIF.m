function y = RK2_EIF(name,v0,para,bool,N,tstep,etime,savedata,reset)  %RK4 method
    if nargin < 9
        reset = false;
        if nargin < 8
            savedata = false;
        end
    end
    y = RK2_IF(name,v0,para,bool,N,tstep,etime,savedata,true,reset);
end
