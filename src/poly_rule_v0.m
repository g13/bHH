function poly_rule_v0(name,picformat,ld,theme,vorigin)
    if nargin < 5
        vorigin = true;
        if nargin < 4
            theme = '';
            if nargin < 3
                ld = true; 
                if nargin < 2
                    picformat = '';
                end
            end
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
    f_I = 5e-7;
    f_E = 5e-7;
    if ld 
        load(['data-v0-',name,theme,'.mat']);
    else
        tstep = 0.02;
        dur = 200;
        vThresh = 0.0;
        dv0 = 0.2;
        vtt = 2.0;
        v0 = dv0:dv0:vtt;
        v0Iterate(tstep,dur,name,f_E,f_I,v0,vThresh,vorigin,theme);
        load(['data-v0-',name,theme,'.mat']);
    end
    N = size(gotEspike,2);
    noEspike = ~gotEspike;
    noIspike = ~gotIspike;
    noSpike = logical(noIspike.*noEspike);
    nValid = sum(noSpike,1);
    power=3; %polynomial power
    tp=12.5;  %observation time
    tp=round(tp/tstep);
    
    ipick = 1:N;
    n = length(ipick);
    h = figure;    
    subplot(2,2,1);
    hold on;
    porders = @(n) n*(n+3)/2;
    CF = zeros(porders(power)-2,n);
    PV = zeros(nv0,n);
    EPSP = PV;
    IPSP = PV;
    SSP = PV;
    residual = PV;
    for i=1:n
        pick = 1:nValid(i);
        E = reshape(VE{ipick(i)}(tp,noSpike(:,i)),[nValid(i),1]);
        I = reshape(VI{ipick(i)}(tp,noSpike(:,i)),[nValid(i),1]);
        S = reshape(VS{ipick(i)}(tp,noSpike(:,i)),[nValid(i),1]);
        EPSP(pick,i) = E;
        IPSP(pick,i) = I;
        SSP(pick,i) = S;
        [CF(:,i),PV(pick,i)]=p_fit11o(E,I,S,power);
        
        plot(PV(pick,i),S,'.','Markersize',15);

        residual(pick,i) = PV(pick,i)-S;
    end
    PV0 =  PV(noSpike);
    minPV = min(PV0(:));
    maxPV = max(PV0(:));
    x=linspace(minPV-0.05*abs(minPV),maxPV+0.05*maxPV,10);
    y=x;
    plot(x,y,'k:','Linewidth',1);
    
    xlabel('Predicted Sum (mV)');
    ylabel('Simulated Sum (mV)'); 
    title(['V_{E} + V_{I} + O(V_{E},V_{I}) ']);
    subplot(2,2,2);
    hold on
    power = 3;
    CF11k = zeros(porders(power)-4,n);
    PV11k = zeros(nv0,n);
    residual11k = PV11k;
    for i=1:n
        pick = 1:nValid(i);
        if ~isempty(pick)
            [CF11k(:,i),PV11k(pick,i)]=p_fit110k0o(EPSP(pick,i),IPSP(pick,i),SSP(pick,i),power);

            p = plot(PV11k(pick,i),SSP(pick,i),'.','Markersize',10);
            plot(EPSP(pick,i)+IPSP(pick,i),SSP(pick,i),'o','Markersize',4,'Color',p.Color);

            residual11k(pick,i) = PV11k(pick,i)-SSP(pick,i);
        end
    end
    title('V_{E} + V_{I} + k_{EI}V_{E}V_{I} + O(V_{E},V_{I})_{>2}');
    
    subplot(2,2,3)
    hold on
    CF2 = zeros(3,n);
    PV2 = zeros(nv0,n);
    residual2 = PV2;
    for i=1:n
        pick = 1:nValid(i);
        if ~isempty(pick)
            [CF2(:,i),PV2(pick,i)]=p_fit110k0(EPSP(pick,i),IPSP(pick,i),SSP(pick,i));

            p = plot(PV2(pick,i),SSP(pick,i),'.','Markersize',10);
            plot(EPSP(pick,i)+IPSP(pick,i),SSP(pick,i),'o','Markersize',4,'Color',p.Color);

            residual2(pick,i) = PV(pick,i)-SSP(pick,i);
        end
    end
    PV20 = PV2(noSpike);
    minV2 = min(PV20);
    maxV2 = max(PV20);
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
        pick = 1:nValid(i);
        color = hsv2rgb(c(i,:));
        x = -EPSP(pick,i).*IPSP(pick,i);
        ypred = EPSP(pick,i)+IPSP(pick,i)-PV2(pick,i);
        y0 = EPSP(pick,i)+IPSP(pick,i)-SSP(pick,i);
        [~,idx] = sort(x);
        plot(x(idx),ypred(idx),'-','Color',color);
        plot(x(idx),y0(idx),'o:','Color',color);
    end
    ylabel('-V_{SC}');
    xlabel('V_{E}V_{I}');
    title(['f_E = ',num2str(f_E,'%g'),' f_I = ',num2str(f_I,'%g')]);
    set(gcf,'color','w');
    if ~isempty(picformat)
        set(gcf,'Renderer','Painters')
        set(gcf,'PaperUnits', 'points','PaperPosition', pPosition);
        fname = ['polyRule_v0','-',name,theme,'.',picformat];
        if strcmp(picformat,'fig')
            saveas(gcf,fname);
        else
            print(gcf,fname,printDriver,dpi);
        end
    end
end
