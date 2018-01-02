picformat = '';
pPosition = [0, 0, 1280, 720];
if ~isempty(picformat)
    if strcmp(picformat,'psc2')
        printDriver = ['-de',picformat];
        picformat = 'eps';
    else
        printDriver = ['-d',picformat];
    end
    dpi = '-r100';
end
tstep = [0.05,0.01,0.005,0.001,0.0005];
parameters;
fCurrent = 0e-4;
n = length(tstep);
dur = 50;
v = zeros(n,1);
%parpool(n);
for i = 1:n
    vVec = trace_plot(tstep(i),dur,fCurrent,'test','',false);
    v(i) = vVec(end);
end
dv = abs(v - v(end));
hold on
figure;
loglog(tstep(1:n-1),dv(1:n-1));
ylabel('error');
xlabel('tstep');

if ~isempty(picformat)
    set(gcf,'Renderer','Painters')
    set(gcf, 'PaperUnits', 'points','PaperPosition', pPosition);
    if strcmp(picformat,'fig')
        saveas(h,['CheckConvergence-',name,'.',picformat]);
    else
        print(h,['CheckConvergence-',name,'.',picformat],printDriver,dpi);
    end
end
%delete(gcp);
