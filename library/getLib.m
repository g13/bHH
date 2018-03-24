addpath(genpath('../src/'));
names = {'RS_exc_Rat','RS_inh_Rat','FS_somato_Rat','IB_somato_GuineaPig','LTS_associa_Cat','TR_somato_Rat'};
dur =   [   300      ,    450     ,      150      ,       700           ,      1000       ,      300      ];
%dtRange0 = [0,2,4,8,12,18,22,24,26,30,60,110,170,230];
dtRange0 = {[0,2,4]};
%tstep = 1/32;
tstep = 1/10;
loadData = false;
singleStored = false;
%singleStored = true;
%loadData = true;
v0 = -0.3:0.1:0.7;
fE = [0.125,0.25] * 3e-5;
fI = fE*2;
rateE = 40;
rateI = 40;
mdur = 1000;
linear0 = true;
theme = 'big';
% for iname = 1:5
for iname = 1:1
    name = names{iname};
    %dtRange = dtRange0{iname};
    %dur = dur{iname}
    dtRange = 0;
    dur = 0;
    % name 
        % 'RS_exc_Rat'
        % 'RS_inh_Rat'
        % 'FS_somato_Rat'
        % 'RS_visual_Ferret'
        % 'FS_visual_Ferret'
        % 'IB_somato_GuineaPig'
        % 'IB_V1_Cat'
        % 'LTS_associa_Cat'
        % 'LTS_somato_Rat' fires at resting potential
        % 'TR_somato_Rat'
    pick = 1;
    model = 'HH';
    picformat = 'png';
    tau_er = 1;
    tau_ed = 3;
    tau_ir = 1;
    tau_id = 5;
    paraNoAdap(model,name,tau_er,tau_ed,tau_ir,tau_id);
    draw = true;
    %draw = false;
    ppp = true;
    %ppp = false;
    npool = 11;
    noAdapV_k4(theme,name,pick,model,picformat,draw,ppp,loadData,npool,v0,fE,fI,singleStored,dur,dtRange,tstep,rateE,rateI,mdur,linear0);
end
