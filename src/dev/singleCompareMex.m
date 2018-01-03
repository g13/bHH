picformat = '';
theme = 'test';
iModel = 0; % 0 for HH 1 for EIF 2 for IF
type = 'RS_exc_Rat';
ith = 1;
switch iModel 
    case 0
        model = 'HH';
    case 1
        model = 'EIF';
    case 2
        model = 'IF';
end
libfile = ['../../library/',type,'-',theme,'-',model,'-',num2str(ith),'th.mat'];
parafile = ['../../library/parameters-',type,'-',theme,'-noAdap-',model,'.mat'];
cutoff = true;
rk = 4; % 2,4
run_t = 500;
ignore_t = 10;
rE = 240; % Hz
rI = 10;
seed = 13;
load(libfile,'dtRange','vRange','sEPSP','sIPSP0','sEPSP0','dur','kVEI','kVIE','kVEE','kVII','fE','fI','nE','ndt','l0');
nt0 = size(sEPSP0,1);
nbt = size(kVEI,1);
ndt = size(kVEI,2);
nv = length(vRange);
ndur = size(sEPSP,1);
edur = l0 - ignore_t;
tstep = dur/(ndur-1);
run_nt = round(run_t/tstep) + 1;
vinit = vRange(2);
[cpuTimeSim,cpuTimeBilinear,outputMat,extIn,inID] = ...
    singleNeuronTestMex(libfile,parafile,ith,rE,rI,iModel,run_t,ignore_t,vinit,rk,cutoff,seed);
[printDriver,dpi,pos0,picformat] = picStd(picformat);
FontSize = 16;
set(0,'DefaultAxesFontSize',FontSize);
set(0,'DefaultTextFontSize',FontSize-2);
fStr = [fE,fI];
t = 0:tstep:run_t;
size(outputMat)
simV = outputMat(:,1);
biV  = outputMat(:,2);
liV  = outputMat(:,3);
gE   = outputMat(:,4);
gI   = outputMat(:,5);
%tin = extIn(:,1);
tin = extIn;
%tID = -extIn(:,2);
tID = -inID;
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
maxV = max([max(simV),max(biV),max(liV)]);
subplot(2,1,1);
plot(t,[simV,biV,liV]);
hold on
vEtar = zeros(Ein,1);
plot(tE,minV*ones(1,Ein),'.r');
for i=1:Ein
    iv = floor(tE(i)/tstep);
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
    iv = floor(tI(i)/tstep);
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
signBiV = sign(sum(biV-simV>0)-sum(biV-simV<0));
errorbar(2,mean(abs(biV-simV))*signBiV,std(abs(biV-simV)));
signLiV = sign(sum(liV-simV>0)-sum(liV-simV<0));
errorbar(1,mean(abs(liV-simV))*signLiV,std(abs(liV-simV)));
xlim([0,3]);
set(gca,'xtick',[1,2]);
plot([0,3],[0,0],':k');
% hold on
% tid = 10;
% vid = find(vRange-vItar(tid)>0,1,'first');
% if ~isempty(vid)
%     plot(1:nt0,sIPSP0(:,Iid(tid),nv),'b');
%     vid = nv;
% else
%     plot(1:nt0,sIPSP0(:,Iid(tid),vid),'r');
%     if vid > 1
%         plot(1:nt0,sIPSP0(:,Iid(tid),vid-1),'b');
%         vid = vid - 1;
%     end
% end
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
% tjd = 9;
% dt = tI(tid) - tI(tjd);
% idt = find(dtRange-dt>0,1,'first');
% if ~isempty(idt)
%     plot(1:nbt,kVII(:,ndt,vid),'b');
% else
%     plot(1:nbt,kVII(:,idt,vid),'r');
%     if idt > 1
%         plot(1:nbt,kVII(:,idt-1,vid),'b');
%     end
% end

% plot(t,gE);
% plot(t,gI);

saveas(gcf,'testMex.fig');
