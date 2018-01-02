function poly_rule(name,picformat,ld)
    if nargin < 3
        ld = true; 
        if nargin < 2
            picformat = '';
        end
    end
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
    if ld 
        load(['data-',name,'.mat']);
    else 
        tstep = 0.01;
        dur = 200;
        f_I = linspace(5e-7,5e-6,5);
        f_E = linspace(3e-8,1e-6,5);
        fInterate(tstep,dur,name,f_E,f_I);
        load(['data-',name,'.mat']);
    end
    power=2; %polynomial power
    tp=25;  %observation time
    tp=round(tp/tstep);
    N = size(VE,3);
    ipick = 1:N;
    n = length(ipick);
    figure(1);    
    subplot(2,1,1);
    hold on;
    porders = @(n) n*(n+3)/2;
    CF = zeros(porders(power),n);
    PV = zeros(nE*nI,n);
    EPSP = PV;
    IPSP = PV;
    SSP = PV;
    residual = PV;
    for i=1:n
        E = reshape(VE(tp,:,ipick(i)),[nE,1]);
        I = reshape(VI(tp,:,ipick(i)),[nI,1]);
        [E,I] = meshgrid(E,I);
        I = reshape(I,[nE*nI,1]);
        E = reshape(E,[nE*nI,1]);
        S = reshape(VS(tp,:,ipick(i)),[nE*nI,1]);
        EPSP(:,i) = E;
        IPSP(:,i) = I;
        SSP(:,i) = S;
        [CF(:,i),PV(:,i)]=p_fit(E,I,S,power);
        
        plot(PV(:,i),S,'.','Markersize',15);

        residual(:,i) = PV(:,i)-S;
    end
    minPV = min(PV(:));
    maxPV = max(PV(:));
    x=linspace(minPV-0.05*abs(minPV),maxPV+0.05*maxPV,10);
    y=x;
    plot(x,y,'k:','Linewidth',1);
    
    xlabel('Predicted Sum (mV)');
    ylabel('Simulated Sum (mV)'); 
    title(['full ',num2str(power),' orders']);
    subplot(2,2,3)
    hold on
    CF2 = zeros(n,1);
    PV2 = zeros(nE*nI,n);
    residual2 = PV2;
    for i=1:n
        [CF2(i),PV2(:,i)]=p_fit2(EPSP(:,i),IPSP(:,i),SSP(:,i));
        
        plot(PV2(:,i),SSP(:,i),'.','Markersize',15);

        residual2(:,i) = PV(:,i)-SSP(:,i);
    end
    minV2 = min(PV2(:));
    maxV2 = max(PV2(:));
    x=linspace(minV2-0.05*abs(minV2),maxV2+0.05*maxV2,10);
    y=x;
    plot(x,y,'k:','Linewidth',1);
    
    xlabel('Predicted Sum (mV)');
    ylabel('Simulated Sum (mV)'); 
    title('V_{E}+V_{I}+k_{EI}*V_{E}V_{I}');
    subplot(2,2,4);
    c= ones(n,3);
    c(:,1) = linspace(0,1,n)';
    c(:,2) = 1;
    c(:,3) = 1;
    hold on
    for i=1:n
        color = hsv2rgb(c(i,:));
        x = -EPSP(:,i).*IPSP(:,i);
        ypred = EPSP(:,i)+IPSP(:,i)-PV2(:,i);
        y0 = EPSP(:,i)+IPSP(:,i)-SSP(:,i);
        [~,idx] = sort(x);
        plot(x(idx),ypred(idx),'-','Color',color);
        plot(x(idx),y0(idx),':','Color',color);
    end
    ylabel('-V_{SC}');
    xlabel('V_{E}V_{I}');
    
    set(gcf,'color','w');
    if ~isempty(picformat)
        set(gcf,'Renderer','Painters')
        set(gcf,'PaperUnits', 'points','PaperPosition', pPosition);
        fname = ['polyRule','-',name,'.',picformat];
        if strcmp(picformat,'fig')
            saveas(gcf,fname);
        else
            print(gcf,fname,printDriver,dpi);
        end
    end
end