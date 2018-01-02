function trace = trace_plot(tstep,f_E,f_I,v0Percent,dur,fCurrent,name,picformat,draw,vonly,returnVonly,silico)
if nargin < 12
    silico = @RK4;
    if nargin < 11
        returnVonly = true;
        if nargin <10
            vonly = false;
            if nargin < 9
                draw = true;
                if nargin < 8
                    picformat = '';
                    if nargin < 7
                        name = 'test';
                        if nargin < 6
                            fCurrent = 1e-4;
                            if nargin < 5
                                dur = 500;
                            end
                        end
                    end
                end
            end    
        end
    end
end
load(['parameters-',name],'para','bool','type','species','posInBrain','n');

pPosition = [0, 0, 1280, 720];
if ~isempty(picformat)
    if strcmp(picformat,'psc2')
        printDriver = ['-de',picformat];
        picformat = 'eps';
    else
        printDriver = ['-d',picformat];
    end
    dpi = '-r300';
end
addpath('./channels/');
para.f_E = f_E;
para.f_I = f_I;
para.fCurrent = fCurrent;
v0 = para.vRest+(para.vT-para.vRest)*v0Percent;
% v0 = para.vLeak;
% v0 = (para.vLeak + para.vT)/2;
% ipick = 1:n;
ipick = 1;
y=silico(name,v0,para,bool,n,tstep,dur);
% y=RK4_0(name,v0,para,bool,n,tstep,dur);
if draw
    t=0:tstep:dur;
    
    for i=1:length(ipick)
        ypick= y(:,:,ipick(i));
        v = ypick(1,:);
        m = ypick(2,:);
        n = ypick(3,:);
        h = ypick(4,:);
        p = ypick(5,:);
        q = ypick(6,:);
        r = ypick(7,:);
        s = ypick(8,:);
        u = ypick(9,:);
        [ge,gi] = synCondTrace(para,t,ipick(i));
        hV = figure;
        hold on      
        ax = plotyy(t,v,t,[ge',gi']);
        xlabel('time (ms)');
        ylabel(ax(1),'voltage (mV)');
        ylabel(ax(2),'conductance (mS/cm^2)');
        legend('v','gE','gI');
        title([type{para.type(ipick(i))},' ',species{para.species(ipick(i))},' ',posInBrain{para.posInBrain(ipick(i))}]);
        if ~isempty(picformat)
            set(hV,'Renderer','Painters')
            set(hV, 'PaperUnits', 'points','PaperPosition', pPosition);
            if strcmp(picformat,'fig')
                saveas(hV,['Volt',num2str(ipick(i)),'-',name,'.',picformat]);
            else
                print(hV,['Volt',num2str(ipick(i)),'-',name,'.',picformat],printDriver,dpi);
            end
        end
        if ~vonly
        hGating = figure; 
        if bool(para.type(ipick(i)),5)
            subplot(2,1,1) 
        end
        hold on;
        box on;
        plot(t,m);
        plot(t,n);
        plot(t,h);
        lg = {'m','n','h'};
        if bool(para.type(ipick(i)),3)
            plot(t,p);
            lg{length(lg)+1} = 'p';
        end
        if bool(para.type(ipick(i)),4)
            plot(t,q);
            plot(t,r);
            lg{length(lg)+1} = 'q';
            lg{length(lg)+1} = 'r';
        end
        legend(lg);
        if bool(para.type(ipick(i)),5)
            subplot(2,2,3) 
            hold on;
            plot(t,s);
            plot(t,u);
            plot(t,u_inf(v,para.vX));
            legend('s','u','u_{\infty}');
            subplot(2,2,4)
            hold on
            plot(t,tau_u(v,para.vX));
            legend('\tau_u');
        end
        
        xlabel('time (ms)');
        ylabel('values of gating variables');
        title([type{para.type(ipick(i))},' ',species{para.species(ipick(i))},' ',posInBrain{para.posInBrain(ipick(i))}]);
        if ~isempty(picformat)
            set(hGating,'Renderer','Painters')
            set(hGating, 'PaperUnits', 'points','PaperPosition', pPosition);
            if strcmp(picformat,'fig')
                saveas(hGating,['gating',num2str(ipick(i)),'-',name,'.',picformat]);
            else
                print(hGating,['gating',num2str(ipick(i)),'-',name,'.',picformat],printDriver,dpi);
            end
        end
        
        hCurrent = figure; hold on
        iLeak = -para.gLeak(ipick(i)).*(v-para.vLeak(ipick(i)));
        iNa = -para.gNa(ipick(i)).*m.^3.*h.*(v-para.vNa);
        iKd = -para.gK(ipick(i)).*n.^4.*(v-para.vK);
        iM = -para.gM(ipick(i)).*p.*(v-para.vK);
        iL = -para.gL(ipick(i)).*q.^2.*r.*(v-para.vCA);
        iT = -para.gT(ipick(i)).*s.^2.*u.*(v-para.vCA);
        
        plot(t,iLeak);
        plot(t,iNa);
        plot(t,iKd);
        plot(t,para.current(t,para.fCurrent/para.S(ipick(i))));
        plot(t,-ge.*(v-para.vE),'r');
        plot(t,-gi.*(v-para.vI),'b');
        lg = {'iLeak','iNa','iKd','iExt','iE','iI'};
        if bool(para.type(ipick(i)),3)
            plot(t,iM);
            lg{length(lg)+1} = 'iM';
        end
        if bool(para.type(ipick(i)),4)
            plot(t,iL);
            lg{length(lg)+1} = 'iL';
        end
        if bool(para.type(ipick(i)),5)
            plot(t,iT);
            lg{length(lg)+1} = 'iT';
        end
        %ylim([]);
        
        legend(lg);
        xlabel('time (ms)');
        ylabel('currents');
        title([type{para.type(ipick(i))},' ',species{para.species(ipick(i))},' ',posInBrain{para.posInBrain(ipick(i))}]);
        
        if ~isempty(picformat)
            set(hCurrent,'Renderer','Painters')
            set(hCurrent,'PaperUnits', 'points','PaperPosition', pPosition);
            if strcmp(picformat,'fig')
                saveas(hCurrent,['Current',num2str(ipick(i)),'-',name,'.',picformat]);
            else
                print(hCurrent,['Current',num2str(ipick(i)),'-',name,'.',picformat],printDriver,dpi);
            end
        end
        end
    end
end
if returnVonly
    trace = y(1,:,ipick);
else
    trace = y;
end
end
