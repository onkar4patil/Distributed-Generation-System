close all
clear all

%constants
rs=25;	 
ls=1*10e-3;



h=1e-6;%step size
j=0;
i=0;
ifilter = 0;
k=0;

rfilter = 10;
cfilter = 1e-4;
lfilter = (1/(120*pi*sqrt(cfilter)))^2;
kfilter = h/(rfilter*cfilter);

kfilter1 = h^2/(lfilter*cfilter);
kfilter2 = h/cfilter;

t=0;
tmax=0.1;

teta1=0;
w1=2*pi*60;
tpwm=0;
Ts=100e-6; %switching period (fs=10kHz)

ias=0;
ibs=0;
ics=0;

Vdc=500;

Vouta = 0;
Voutb = 0;
Voutc = 0;

while t<=tmax,

    t = t + h; % simulation time

    teta1=teta1+h*w1; 
    if teta1 >= 2*pi,
        teta1=teta1-2*pi;
    end
    
    vas_ref = 25*cos(teta1);
    vbs_ref = 25*cos(teta1-2*pi/3);
    vcs_ref = 25*cos(teta1-4*pi/3);
    
    %begin of the pwm loop simulation (interruption)
    if t>tpwm,
        
        tpwm=tpwm+Ts;
        
        cont=0;
        
        tal1=(vas_ref/Vdc+0.5)*Ts;
        tal2=(vbs_ref/Vdc+0.5)*Ts;        
        tal3=(vcs_ref/Vdc+0.5)*Ts;
        
    end
    %end of the pwm loop simulation

    cont=cont+h;

    %begin simulation of the converter
    if (cont<tal1)
        q1=1;
    else
        q1=0;
    end    
    if (cont<tal2)
        q2=1;
    else
        q2=0;
    end     
    if (cont<tal3)
        q3=1;
    else
        q3=0;
    end     

    %Definition of the pole voltages
    v10 = (2*q1-1)*(Vdc/2);
    v20 = (2*q2-1)*(Vdc/2);
    v30 = (2*q3-1)*(Vdc/2);
    
    vn0 = (1/3)*(v10+v20+v30);
    
    %Phase voltages applied to the stator of the machine
    vas = v10 - vn0;
    vbs = v20 - vn0;
    vcs = v30 - vn0;

%    Vouta = Vouta + h^2/(lfilter*cfilter+h^2*lfilter*cfilter)*(v10-Vouta);
%    Voutb = Voutb + h^2/(lfilter*cfilter+h^2*lfilter*cfilter)*(v20-Voutb);
%    Voutc = Voutc + h^2/(lfilter*cfilter+h^2*lfilter*cfilter)*(v30-Voutc);
    
    
%    Vouta = Vouta + kfilter1*(v10-Vouta) - kfilter2*Vouta;
%    Voutb = Voutb + kfilter1*(v20-Voutb) - kfilter2*Voutb;
%    Voutc = Voutc + kfilter1*(v30-Voutc) - kfilter2*Voutc;
    
    
    Vouta = Vouta + kfilter*(v10-Vouta);
    Voutb = Voutb + kfilter*(v20-Voutb);
    Voutc = Voutc + kfilter*(v30-Voutc);

    ias = ias + (Vouta - rs*ias)*h/ls;
    ibs = ibs + (Voutb - rs*ibs)*h/ls;
    ics = ics + (Voutc - rs*ics)*h/ls;
    
    P = Vouta*ias+Voutb*ibs+Voutc*ics;

    % storage the variables
    
    j = j + 1;
    
    time(j)=t;
    
    voltage_as(j)=vas;
    voltage_bs(j)=vbs;
    voltage_cs(j)=vcs;    

    voltage_as_ref(j)=vas_ref;
    voltage_bs_ref(j)=vbs_ref;
    voltage_cs_ref(j)=vcs_ref;       
    
    current_as(j)=ias;
    current_bs(j)=ibs;
    current_cs(j)=ics;
    
    pole_v10(j)=v10;
    pole_v20(j)=v20;
    pole_v30(j)=v30;
    
    voutas(j)=Vouta;
    voutbs(j)=Voutb;
    voutcs(j)=Voutc;
    
    Power(j) = P;

  
    
end


figure(1),plot(time,voltage_as_ref,time,voltage_bs_ref,time,voltage_cs_ref)

figure(2),plot(time,pole_v10,time,pole_v20,time,pole_v30)

figure(3),plot(time,voltage_as,time,voltage_bs,time,voltage_cs)

figure(4),plot(time,current_as,time,current_bs,time,current_cs)

figure(5),plot(time,voutas,time,voutbs,time,voutcs)

figure(6),plot(time,Power)





