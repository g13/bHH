format long
picformat = '';
theme = 'test';
iModel = 0; % 0 for HH 1 for EIF 2 for IF
type = 'RS_exc_Rat';
ith = 1;
tref = 7;
switch iModel 
    case 0
        model = 'HH';
    case 1
        model = 'EIF';
    case 2
        model = 'IF';
end
libfile = ['../../library/',type,'-',theme,'-',model,'-',num2str(ith),'th.mat'];
parafile = ['../../library/parameters-',type,'-noAdap-',model,'.mat'];
cutoff = false;
rk = 4; % 2,4
run_t = 1000;
ignore_t = 10;
% rE = [15,20,25,30,35,40,45]; % Hz
% rI = [30,35,40,45,50,55,60];
rE = [70];
rI = [40];
rEl = length(rE);
assert(rEl==length(rI));
seed = 13;
load(libfile,'dtRange','vRange','sEPSP','sIPSP0','sEPSP0','dur','kV','kEI','kIE','kEE','kII','fE','fI','nE','ndt','l0');
nt0 = size(sEPSP0,1)
nt = size(kEI,1)
ndt = size(kEI,2)
nv = length(vRange)
ndur = size(sEPSP,1)
edur = l0 - ignore_t;
tstep = dur/(ndur-1)
run_nt = round(run_t/tstep) + 1;
vinit = vRange(2);
afterSpikeBehavior = 2;
spikeShape = true;
kVStyle = true;
[cpuTimeAndTsp,outputMats,extIns,inIDs,tspCell] = ...
    singleNeuronMex(libfile,parafile,ith,rE,rI,iModel,run_t,ignore_t,vinit,rk,cutoff,seed,tref,afterSpikeBehavior,spikeShape,kVStyle);
[printDriver,dpi,pos0,picformat] = picStd(picformat);
FontSize = 16;
set(0,'DefaultAxesFontSize',FontSize);
set(0,'DefaultTextFontSize',FontSize-2);
fStr = [fE,fI];
t = 0:tstep:run_t;
cpuTime = cpuTimeAndTsp(1:3,:);
tspSize = cpuTimeAndTsp(4:6,:);
err = zeros(rEl,2,2);
%
for j=1:rEl
    outputMat = outputMats(:,:,j);
    extIn = extIns{j};
    tID = inIDs{j};
    tsp_sim = tspCell{j};
    tsp_bi = tspCell{j+1};
    tsp_li = tspCell{j+2};
    simV = outputMat(:,1);
    biV  = outputMat(:,2);
    liV  = outputMat(:,3);
    gE   = outputMat(:,4);
    gI   = outputMat(:,5);
    if iModel == 0
        m = outputMat(:,6);
        n = outputMat(:,7);
        h = outputMat(:,8);
    end
    tin = extIn;
    tID = tID+1;
    pE = (tID <= nE);
    pI = (tID > nE);
    tE = tin(pE);
    tI = tin(pI);
    Eid = tID(pE);
    Iid = tID(pI)-nE;
    Ein = length(tE);
    Iin = length(tI);
    figure;
    textFontSize = 8;
    minV = min([min(simV),min(biV),min(liV)]);
    maxV = min([-55,max([max(simV),max(biV),max(liV)])]);
    subplot(2,1,1);
    hold on
    plot(t,[simV,biV,liV]);
    plot(t,gE,':r');
    plot(t,gI,':b');
    %plot(t,[simV]);
    xlabel('t ms');
    ylabel('mem mV');
    hold on
    vEtar = zeros(Ein,1);
    plot(tE,minV*ones(1,Ein),'.r');
    for i=1:Ein
        iv = floor(tE(i)/tstep)+1;
        jv = iv + 1;
        vEtar = simV(iv) + mod(tE(i),tstep)/tstep * (simV(jv)-simV(iv));
        plot([tE(i),tE(i)],[minV,vEtar],':r');
        text(tE(i),minV+(vEtar-minV)*0.3,num2str(Eid(i)),'Color','r','FontSize',textFontSize);
        text(tE(i),minV+(vEtar-minV)*0.1,num2str(i),'Color','r','FontSize',textFontSize);
     
        iv = floor((tE(i)+edur)/tstep);
        if iv+1 < run_nt
            jv = iv + 1;
            vEtar = simV(iv) + mod(tE(i),tstep)/tstep * (simV(jv)-simV(iv));
            plot([tE(i),tE(i)]+edur,[vEtar,maxV],':r');
            text(tE(i)+edur,maxV-(maxV-vEtar)*0.3,num2str(Eid(i)),'Color','r','FontSize',textFontSize);
            text(tE(i)+edur,maxV-(maxV-vEtar)*0.1,num2str(i),'Color','r','FontSize',textFontSize);
        end
    end
    vItar = zeros(Iin,1);
    plot(tI,minV*ones(1,Iin),'.b');
    for i=1:Iin
        iv = floor(tI(i)/tstep)+1;
        jv = iv + 1;
        vItar(i) = simV(iv) + mod(tI(i),tstep)/tstep * (simV(jv)-simV(iv));
        plot([tI(i),tI(i)],[minV,vItar(i)],':b');
        text(tI(i),minV+(vItar(i)-minV)*0.3,num2str(Iid(i)),'Color','b','FontSize',textFontSize);
        text(tI(i),minV+(vItar(i)-minV)*0.1,num2str(i),'Color','b','FontSize',textFontSize);
            
        iv = floor((tI(i)+edur)/tstep);
        if iv+1 < run_nt
            jv = iv + 1;
            vtar = simV(iv) + mod(tI(i),tstep)/tstep * (simV(jv)-simV(iv));
            plot([tI(i),tI(i)]+edur,[vtar,maxV],':b');
            text(tI(i)+edur,maxV-(maxV-vtar)*0.3,num2str(Iid(i)),'Color','b','FontSize',textFontSize);
            text(tI(i)+edur,maxV-(maxV-vtar)*0.1,num2str(i),'Color','b','FontSize',textFontSize);
        end
    end
    legend({'sim','bi','li'});
    ylim([minV,maxV]);
    subplot(2,2,3);
    hold on;
    signLiV = sign(sum(liV-simV>0)-sum(liV-simV<0));
    err(j,1,1) = mean(abs(liV-simV))*signLiV;
    err(j,2,1) = std(abs(liV-simV));
    errorbar(1,err(j,1,1),err(j,2,1));
    signBiV = sign(sum(biV-simV>0)-sum(biV-simV<0));
    err(j,1,2) = mean(abs(biV-simV))*signBiV;
    err(j,2,2) = std(abs(biV-simV));
    errorbar(2,err(j,1,2),err(j,2,2));
    xlim([0,3]);
    set(gca,'xtick',[1,2]);
    plot([0,3],[0,0],':k');
    set(gca,'xtick',[1,2]);
    set(gca,'xticklabel',{'linear','bilinear'});
    ylabel('error per timestep mV');
    subplot(2,2,4);
    dliV = liV - simV;
    plot(t,dliV);
    hold on
    dbiV = biV - simV;
    plot(t,dbiV);
    legend({'\Delta(sim-li)','\Delta(sim-bi)'});
    minV = min(min(dbiV),min(dliV));
    maxV = max(max(dbiV),max(dliV));
    vEtar = zeros(Ein,1);
    plot(tE,minV*ones(1,Ein),'.r');
    targetV = dliV;
    for i=1:Ein
        iv = floor(tE(i)/tstep);
        jv = iv + 1;
        vEtar = targetV(iv) + mod(tE(i),tstep)/tstep * (targetV(jv)-targetV(iv));
        plot([tE(i),tE(i)],[minV,vEtar],':r');
        text(tE(i),minV+(vEtar-minV)*0.3,num2str(Eid(i)),'Color','r','FontSize',textFontSize);
        text(tE(i),minV+(vEtar-minV)*0.1,num2str(i),'Color','r','FontSize',textFontSize);
     
        iv = floor((tE(i)+edur)/tstep);
        if iv+1 < run_nt
            jv = iv + 1;
            vEtar = targetV(iv) + mod(tE(i),tstep)/tstep * (targetV(jv)-targetV(iv));
            plot([tE(i),tE(i)]+edur,[vEtar,maxV],':r');
            text(tE(i)+edur,maxV-(maxV-vEtar)*0.3,num2str(Eid(i)),'Color','r','FontSize',textFontSize);
            text(tE(i)+edur,maxV-(maxV-vEtar)*0.1,num2str(i),'Color','r','FontSize',textFontSize);
        end
    end
    vItar = zeros(Iin,1);
    plot(tI,minV*ones(1,Iin),'.b');
    for i=1:Iin
        iv = floor(tI(i)/tstep);
        jv = iv + 1;
        vItar(i) = targetV(iv) + mod(tI(i),tstep)/tstep * (targetV(jv)-targetV(iv));
        plot([tI(i),tI(i)],[minV,vItar(i)],':b');
        text(tI(i),minV+(vItar(i)-minV)*0.3,num2str(Iid(i)),'Color','b','FontSize',textFontSize);
        text(tI(i),minV+(vItar(i)-minV)*0.1,num2str(i),'Color','b','FontSize',textFontSize);
            
        iv = floor((tI(i)+edur)/tstep);
        if iv+1 < run_nt
            jv = iv + 1;
            vtar = targetV(iv) + mod(tI(i),tstep)/tstep * (targetV(jv)-targetV(iv));
            plot([tI(i),tI(i)]+edur,[vtar,maxV],':b');
            text(tI(i)+edur,maxV-(maxV-vtar)*0.3,num2str(Iid(i)),'Color','b','FontSize',textFontSize);
            text(tI(i)+edur,maxV-(maxV-vtar)*0.1,num2str(i),'Color','b','FontSize',textFontSize);
        end
    end
    saveas(gcf,['testMex-',num2str(j),'-E',num2str(rE(j)),'-I',num2str(rI(j)),'.fig']);
end
figure;
subplot(1,2,1)
plot(rE+rI,cpuTime');
legend({'rk2','bilinear','linear'});
xlabel('E+I input rate');
ylabel('cpu time cost');
subplot(1,2,2)
hold on
errorbar(rE+rI,err(:,1,1),err(:,2,1),'b');
errorbar((rE+rI)+0.2,err(:,1,2),err(:,2,2),'r');
xlabel('E+I input rate');
ylabel('error per timestep vs RK2');
legend({'linear','bilinear'});
saveas(gcf,['err_rate','.fig']);
