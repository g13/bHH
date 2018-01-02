function paraNoAdap(model,name,theme,tau_er,tau_ed,tau_ir,tau_id,draw)
    if nargin < 8
        draw = false;
    end
%%
%    mh,  n, p,qr,su
%   iNa,iKd,iM,iL,iT
% RS  1   1  1  0  0
% FS  1   1  1  0  0
% IB  1   1  1  1  0
% LTS 1   1  1  0  1
% TR  1   1  0  0  0.5
bool = [1,1,1,0,0;
        1,1,1,0,0;
        1,1,1,1,0;
        1,1,1,0,1;
        1,1,0,0,1];
% TR's gT is modelled by a nonlinear constant field equations. Destexhe et al., Hille 1992.
type = {'RS','FS','IB','LTS','TR'};

species = {'rat','guinea pig','cat','ferret'};
% whole-cell patch for rat, sharp electrode for others. 
% cat in vivo. others in vitro.
                %   1         ,         2        ,       3       ,         4       ,        5          6    %
posInBrain = {'somato. cortex','somato. thalamus','visual cortex','associa. cortex','frontal cortex','V1'};

%           L=d,    vleak,      vT,     Vx,     gleak,     gNa,    gKd,    gM,     gL,     gT,     tau_max,   TYPE,    species,    posInBrain
dataSet = [... %   RS exc rat somato. cortex   %
           61.4,    -70.3,   -59.72,   inf,     2.73,       39,    6.0,    0.20,   inf,    inf,    1445.0,     1,      1,          1,           1;
           61.4,    -70.3,   -56.16,   inf,     2.43,       56,    6.0,   0.075,   inf,    inf,     608.0,     1,      1,          1,           1;
           61.4,    -70.3,   -65.23,   inf,     1.33,       60,    5.1,   0.087,   inf,    inf,    2269.0,     1,      1,          1,           1;
           61.4,    -70.3,   -55.43,   inf,     1.94,       52,    3.7,    0.15,   inf,    inf,     653.5,     1,      1,          1,           1;
           61.4,    -70.3,   -62.86,   inf,     2.41,       36,    3.1,   0.088,   inf,    inf,    1476.0,     1,      1,          1,           1;
           61.4,    -70.3,   -62.14,   inf,     1.43,       50,    6.0,   0.097,   inf,    inf,     932.2,     1,      1,          1,           1;
           61.4,    -70.3,   -66.51,   inf,     1.42,       60,    5.5,    0.20,   inf,    inf,    1340.0,     1,      1,          1,           1;
           61.4,    -70.3,   -62.87,   inf,     1.04,       58,    5.9,    0.10,   inf,    inf,     959.0,     1,      1,          1,           1;
           61.4,    -70.3,   -59.54,   inf,     0.98,       58,    6.0,   0.082,   inf,    inf,    1351.0,     1,      1,          1,           1;
           61.4,    -70.3,   -62.96,   inf,     0.91,       45,    2.0,    0.10,   inf,    inf,     583.3,     1,      1,          1,           1;
           61.4,    -70.3,   -58.67,   inf,     2.23,       59,    3.1,    0.16,   inf,    inf,     686.4,     1,      1,          1,           1;
           61.4,    -70.3,   -63.94,   inf,     2.54,       42,    3.9,    0.20,   inf,    inf,     610.9,     1,      1,          1,           1;
           61.4,    -70.3,   -63.62,   inf,     3.34,       30,    6.0,    0.10,   inf,    inf,    1691.0,     1,      1,          1,           1;
           61.4,    -70.3,   -67.42,   inf,     3.06,       40,    5.7,    0.20,   inf,    inf,    2928.0,     1,      1,          1,           1;
            ... %   RS Inh rat somato. cortex, vLeak unknown  %                                                                     
           61.8,    -76.2,   -61.89,   inf,     1.71,       27,    2.6,   0.038,   inf,    inf,    1327.0,     1,      1,          1,           0;
           61.8,    -76.2,   -60.01,   inf,     0.96,       36,    5.5,   0.073,   inf,    inf,    1996.0,     1,      1,          1,           0;
           61.8,    -76.2,   -74.67,   inf,     0.70,       18,    4.2,    0.17,   inf,    inf,    1541.0,     1,      1,          1,           0;
           61.8,    -76.2,   -71.66,   inf,     5.69,       14,    5.2,    0.20,   inf,    inf,    1490.0,     1,      1,          1,           0;
           61.8,    -76.2,   -66.54,   inf,     1.12,       28,    6.0,    0.09,   inf,    inf,    2646.0,     1,      1,          1,           0;
           61.8,    -76.2,   -59.29,   inf,     1.83,       40,    3.3,   0.017,   inf,    inf,    1594.0,     1,      1,          1,           0;
           61.8,    -76.2,   -62.51,   inf,     1.46,       10,    7.0,   0.035,   inf,    inf,    1349.0,     1,      1,          1,           0;
           61.8,    -76.2,   -62.37,   inf,     2.19,       40,    2.9,   0.044,   inf,    inf,    2997.0,     1,      1,          1,           0;
           61.8,    -76.2,   -67.85,   inf,     1.60,       10,    2.1,   0.098,   inf,    inf,     934.4,     1,      1,          1,           0;
            ... %   FS rat somato. cortex   %                                                                                       
           56.9,    -70.4,   -56.29,   inf,    12.35,       31,    7.0,   0.049,   inf,    inf,     505.5,     2,      1,          1,           0;
           56.9,    -70.4,   -64.13,   inf,    10.00,       32,    6.0,    0.09,   inf,    inf,     508.3,     2,      1,          1,           0;
           56.9,    -70.4,   -62.15,   inf,    14.81,       38,    5.2,   0.097,   inf,    inf,     500.0,     2,      1,          1,           0;
           56.9,    -70.4,   -67.65,   inf,    29.06,       51,    5.7,     0.1,   inf,    inf,     500.3,     2,      1,          1,           0;
           56.9,    -70.4,   -58.71,   inf,    12.68,       60,    4.5,     0.1,   inf,    inf,     500.5,     2,      1,          1,           0;
           56.9,    -70.4,   -65.42,   inf,     5.84,       40,    5.4,   0.038,   inf,    inf,     664.2,     2,      1,          1,           0;
           56.9,    -70.4,   -61.47,   inf,     3.86,       58,    6.6,    0.05,   inf,    inf,    1056.0,     2,      1,          1,           0;
           56.9,    -70.4,   -58.20,   inf,     2.66,       44,    6.0,     0.1,   inf,    inf,     503.9,     2,      1,          1,           0;
           56.9,    -70.4,   -63.44,   inf,     5.64,       43,    4.4,     0.1,   inf,    inf,     505.7,     2,      1,          1,           0;
           56.9,    -70.4,   -57.90,   inf,     3.865,      58,    3.9,  0.0787,   inf,    inf,     502.0,     2,      1,          1,           0;
           56.9,    -70.4,   -63.85,   inf,     4.50,       43,    4.7,    0.05,   inf,    inf,    1723.0,     2,      1,          1,           0;
           56.9,    -70.4,   -58.79,   inf,     0.46,       60,    3.9,   0.021,   inf,    inf,    1454.0,     2,      1,          1,           0;
           56.9,    -70.4,   -67.17,   inf,     1.27,       32,    2.2,   0.071,   inf,    inf,     596.6,     2,      1,          1,           0;
           56.9,    -70.4,   -60.49,   inf,     7.95,       55,    6.0,   0.039,   inf,    inf,    2023.0,     2,      1,          1,           0;
... %*** vT for ferret, guinea-pig and cat associa. cortex is an initial guess, not provided in the paper %                         
            ... %   RS ferret visual cortex   %                                                                                     
           96.0,    -70.0,   -61.0,    inf,    28.95,       50,    5.0,    0.07,   inf,    inf,    4000.0,     1,      4,          3,           1;
            ... %   FS ferret visual cortex %                                                                                       
           67.0,    -70.0,   -61.5,    inf,    43.43,       50,   10.0,    0.07,   inf,    inf,    4000.0,     2,      4,          3,           0;
            ... %   IB guinea-pig somato. cortex    %                                                                               
           96.0,    -70.0,   -61.5,    inf,     2.90,       50,    5.0,    0.03,   0.1,    inf,    4000.0,     3,      2,          1,           -1;
           96.0,    -70.0,   -61.5,    inf,     2.90,       50,    5.0,    0.03,   0.2,    inf,    4000.0,     3,      2,          1,           -1;
            ... %   IB cat V1   %                                                                                                   
           96.0,    -75.0,   -58.0,    inf,    28.95,       50,    4.2,   0.042,  0.12,    inf,    1000.0,     3,      3,          6,           -1;
            ... %   LTS cat associa. cortex %                                                                                       
           96.0,    -75.0,   -58.0,      2,     2.90,       50,    5.0,    0.03,   inf,    0.4,    4000.0,     4,      3,          4,           -1;
            ... %   LTS rat somato. cortex  %                                                                                       
           89.2,    -50.0,   -50.0,     -7,     4.75,       50,    4.0,   0.028,   inf,    0.4,    4000.0,     4,      1,          1,           -1;
            ... %   TR rat ventrobasal nucleus Destexhe et al.,  %                                                                  
            ... %  modified from L*d=100*76.6 C=0.878; E_K = -100(-90 for others, but let's use -90 first);                         
           87.5,   -69.85,   -52.0,     -1,    10.39,    113.9,  113.9,     inf,   inf,    0.4,       inf,     5,      1,          2,           -1];
%           L=d,    vleak,      vT,     Vx,     gleak,     gNa,    gKd,    gM,     gL,     gT,     tau_max,   TYPE,    species,    posInBrain
%           1         2         3       4         5         6       7       8       9       10      11          12      13          14
% noAdap = true;
load(['vRest-noAdap-',model,'.mat']);
dataSet = [dataSet, vRest];
N = size(dataSet,1);
para = struct();
    switch name 
        case 'RS_exc_Rat'
            select = find(abs(dataSet(:,1) - 61.4) < 1e-15 & dataSet(:,12) == 1 & dataSet(:,13) == 1 & dataSet(:,14) == 1);
     
        case 'RS_inh_Rat'
            select = find(abs(dataSet(:,1) - 61.8) < 1e-15 & dataSet(:,12) == 1 & dataSet(:,13) == 1 & dataSet(:,14) == 1);
 
        case 'FS_somato_Rat'
            select = find(abs(dataSet(:,1) - 56.9) < 1e-15 & dataSet(:,12) == 2 & dataSet(:,13) == 1 & dataSet(:,14) == 1);

        case 'RS_visual_Ferret'
            select = find(abs(dataSet(:,1) - 96) < 1e-15 & dataSet(:,12) == 1 & dataSet(:,13) == 4 & dataSet(:,14) == 3);
   
        case 'FS_visual_Ferret'
            select = find(abs(dataSet(:,1) - 67) < 1e-15 & dataSet(:,12) == 2 & dataSet(:,13) == 4 & dataSet(:,14) == 3);

        case 'IB_somato_GuineaPig'
            select = find(abs(dataSet(:,1) - 96) < 1e-15 & dataSet(:,12) == 3 & dataSet(:,13) == 2 & dataSet(:,14) == 1);

        case 'IB_V1_Cat'
            select = find(abs(dataSet(:,1) - 96) < 1e-15 & dataSet(:,12) == 3 & dataSet(:,13) == 3 & dataSet(:,14) == 6);

        case 'LTS_associa_Cat'
            select = find(abs(dataSet(:,1) - 96) < 1e-15 & dataSet(:,12) == 4 & dataSet(:,13) == 3 & dataSet(:,14) == 4);

        case 'LTS_somato_Rat'
            select = find(abs(dataSet(:,1) - 89.2) < 1e-15 & dataSet(:,12) == 4 & dataSet(:,13) == 1 & dataSet(:,14) == 1);

        case 'TR_somato_Rat'
            select = find(abs(dataSet(:,1) - 87.5) < 1e-15 & dataSet(:,12) == 5 & dataSet(:,13) == 1 & dataSet(:,14) == 2);
    end
fname = [name,'-',theme,'-noAdap-',model];
para.select = select;

para.S = pi*(dataSet(select,1)*1e-4).^2; % verified by Fig.2 in ref minimal HH, only the side area of the cylinder are included.
n = length(select);
para.gNa = dataSet(select,6); 
para.gK  = dataSet(select,7);
para.gLeak = dataSet(select,5)./para.S*1e-6;
para.gL = dataSet(select,9);
%para.gM = dataSet(select,8);
para.gM = zeros(size(para.gLeak));
para.ei = dataSet(select,15);
para.DeltaT = 0.4375*ones(n,1);
para.gT = dataSet(select,10);
para.tE = [0];
para.f_E = 0e-4;
para.nEspike = length(para.tE);
para.tI = [0];
para.nIspike = length(para.tI);
para.f_I = 0e-4;
para.vRest = dataSet(select,16);
para.vX = dataSet(select,4);
para.vLeak = dataSet(select,2);
para.vNa = 50;
para.vK = -90;
para.vCA = 120;
para.vE = 0;
para.vI = -80;
para.vT = dataSet(select,3);
para.tau_max = dataSet(select,11);
para.type = dataSet(select,12);
para.species = dataSet(select,13);
para.posInBrain = dataSet(select,14);

para.tau_e=tau_ed;
para.tau_er=tau_er;
% this form have g(tmax) = f_E;
% para.gConE=1/((para.tau_er/para.tau_e)^(para.tau_er/(para.tau_e-para.tau_er))-(para.tau_er/para.tau_e)^(para.tau_e/(para.tau_e-para.tau_er)));
para.gConE = 1/(para.tau_e-para.tau_er);
para.tau_i=tau_id;
para.tau_ir=tau_ir;
% para.gConI=1/((para.tau_ir/para.tau_i)^(para.tau_ir/(para.tau_i-para.tau_ir))-(para.tau_ir/para.tau_i)^(para.tau_i/(para.tau_i-para.tau_ir)));
para.gConI = 1/(para.tau_i-para.tau_ir);
%%% When checking convergence rate of the algorithm, ct0 < start time,
% ct1>end time
ct0 = 200;
ct1 = 700;
pwidth = 20;
step = true;
% step = false;
pulse = false;
% pulse = true;
% para.current = @(t,fCurrent) fCurrent./para.S;
%%% step current
if step
    para.current = @(t,fCurrent) (1+sign(t-ct0))/2.*fCurrent + (1-sign(t-ct0))/2.*0.0 - (1+sign(t-ct1))/2.*fCurrent;
end
%%% pulse current
if pulse
    para.current = @(t,fCurrent)(1-sign(t-ct1))/2.*(1+sign(t-ct0))/2.*(1-mod(floor((t-50)/pwidth),2))*fCurrent;
end
save(['parameters-',fname,'.mat'],'para','bool','type','species','posInBrain','n');
n
    if draw
        vonly = false;
        returnVonly = false;
        cplot = false;
        % cplot = true;
        draw0 = true;
        % draw0 = false;
        f_E = 0e-6;
        para.tE = 0;
        f_I = 0e-6;
        para.tI = 0;
        v0Percent = 0.0;
        fCurrent = 1e-4;
        trace_plot(0.05,f_E,f_I,v0Percent,1000,fCurrent,name,'',draw0,vonly,returnVonly,@RK4,cplot,true);
    end
end
