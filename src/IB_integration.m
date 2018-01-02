% scan input strengthes and record the corresponding EPSP, IPSP and SSP
% save the data for processing in polynomial_rule.m

clear;
close all;
clc;
 
initialization();

I_current=0;

x0=fcurrent(vrest);
dur=100;

I_Count=1;
E_Count=1;
S_Count=1;

f_E=0;
for f_I=linspace(5e-6/S,5e-5/S,10)
    x=RK4_IB(@voltagev,@gatingm,@gatingn,@gatingh,@gatingp,@gatingq,@gatingr,x0,dur);
    IPSP(I_Count,:)=x(1,:);
    I_Count=I_Count+1;
end

f_I=0;
for f_E=linspace(3e-7/S,1e-5/S,10)
    x=RK4_IB(@voltagev,@gatingm,@gatingn,@gatingh,@gatingp,@gatingq,@gatingr,x0,dur);
    EPSP(E_Count,:)=x(1,:);
    E_Count=E_Count+1;
end
 
for f_E=linspace(3e-7/S,1e-5/S,10)    
    for f_I=linspace(5e-6/S,5e-5/S,10)
        
        x=RK4_IB(@voltagev,@gatingm,@gatingn,@gatingh,@gatingp,@gatingq,@gatingr,x0,dur);
        
        SSP(S_Count,:)=x(1,:);
        S_Count=S_Count+1;
        disp(S_Count);
     end
end

VI=IPSP(:,1:round(1/tstep):end)-vrest;
VE=EPSP(:,1:round(1/tstep):end)-vrest;
VS=SSP(:,1:round(1/tstep):end)-vrest;

save IB_integration vrest tstep VI VE VS


