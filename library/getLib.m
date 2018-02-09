addpath(genpath('../src/'));
names = {'RS_exc_Rat','RS_inh_Rat','FS_somato_Rat','IB_somato_GuineaPig','LTS_associa_Cat','TR_somato_Rat'};
dur =   [   300      ,    450     ,      150      ,       700           ,      1000       ,      300      ];
dtRange0 = [0,2,4,8,12,18,22,24,26,30,60,110,170,230];
dur0 = 300;
dtRatio = dtRange0./dur0;
tstep = 1/32;
for iname = 1:length(names) 
    name = names{iname};
    dtRange = dtRatio * dur(iname);
    dtRange = round(dtRange/tstep)*tstep;
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
    theme = 'test';
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
    npool = 16;
    loadData = false;
    singleStored = false;
    %singleStored = true;
    %loadData = true;
    v0 = -0.3:0.15:1.2;
    %v0 = [-0.4,0,0.4,1.2];
    %v0 = [-0.3,0,0.5];
    fE = linspace(0.125,0.25,4) * 2e-5;
    fE = fE(1:2);
    %fI = (0.25:0.25:1.0) * 1e-5;
    fI = fE(1:2)*0.5;
    avoidSpike = true;
    %v0 = 0.0:0.1:0.2;
    %noAdapV(theme,name,pick,model,picformat,draw,ppp,loadData,npool,v0,fE,fI,avoidSpike);
    noAdapV_k4(theme,name,pick,model,picformat,draw,ppp,loadData,npool,v0,fE,fI,singleStored,dur(iname),dtRange,tstep);
end
