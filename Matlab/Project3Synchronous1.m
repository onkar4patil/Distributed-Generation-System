%% ECE 610 / Project 3

%% Update 11/13
%% Finished simulation for synchronous machine


%% Questions for Dos Santos 
%%
%% 1) Are our assumptions okay?
%%     a) used motor parameters from table 5.10-2
%%     b) constant torque means constant rotor speed
%%     c) number of poles = 2 for simplification

%% 2) what do we expect from synchronous generator output alone?
%%     a) We have a output of 50MVA. Current lags voltage by ~= pi/4
%%     b) Should it be in phase? or does filter control that on load side
%%     c) What would we expect from the output of the generator

%% Define constants for the synchronous machine described on p176

close all
clear all

index = 1;          %% Define index for vectors
j = sqrt(-1);       %% Define j to be imaginary
Vs = 38.35/sqrt(2);    %% Initial value of V
rs = 0.00243;       %% Define 

Xls = 0.1538;
Xq = 1.457;
Xmq = Xq - Xls;
Xd = 1.457;
Xmd = Xd - Xls;
rkq1 = 0.00144;
rfd = 0.00075;
Xlkq1 = 0.6578;
Xkq1 = Xlkq1 + Xmq;
Xlfd = 0.1145;
Xfd = Xlfd + Xmd;
rkq2 = 0.00681;
rkd = 0.01080;
Xlkq2 = 0.07062;
Xkq2 = Xlkq2 + Xmq;
Xlkd = 0.06577;
Xkd = Xlkd + Xmd;

p = 2;

theta = 0;          %% Initial theta of 0
omega_e = 120*pi;      %% Electrical Speed in rads/sec
omega_r = 120*pi;


%% We should also define some more values for our iterative solution

t = 0; %% Start taking data from t = 0 seconds
h = 10e-6; % step size is 1 microsecond
tmax = 0.1; % Stop at t = 0.1 seconds
Ks(1,:) = [cos(theta), cos(theta-2*pi/3), cos(theta+2*pi/3)];
Ks(2,:) = [sin(theta), sin(theta-2*pi/3), sin(theta+2*pi/3)];
Ks(3,:) = [1/2, 1/2, 1/2];
Ks = (2/3)*Ks;              %% Reference frame Matrix


%% Now to define initial values for Vabc, Iabc, Vqdo, and Iqdo. 

Vabcs(1,1) = sqrt(2)*Vs*cos(omega_e*t);
Vabcs(2,1) = sqrt(2)*Vs*cos(omega_e*t-2*pi/3);
Vabcs(3,1) = sqrt(2)*Vs*cos(omega_e*t+2*pi/3);

Vqdos = Ks * Vabcs;
Vkq1 = 0;
Vkq2 = 0;
exfd = 1;
Vkd = 0;

Vref(1,1) = Vqdos(1,1);        
Vref(2,1) = Vqdos(2,1);
Vref(3,1) = Vqdos(3,1);
Vref(4,1) = Vkq1;
Vref(5,1) = Vkq2;
Vref(6,1) = exfd;
Vref(7,1) = Vkd;

%% Need to build this matrix to find the currents

Zref(1,:) = [rs+(p/omega_e)*Xq, (omega_r/omega_e)*Xd, 0, (p/omega_e)*Xmq, (p/omega_e)*Xmq, (omega_r/omega_e)*Xmd, (omega_r/omega_e)*Xmd];
Zref(2,:) = [-(omega_r/omega_e)*Xq, rs+(p/omega_e)*Xd, 0, -(omega_r/omega_e)*Xmq, -(omega_r/omega_e)*Xmq, (p/omega_e)*Xmd, (p/omega_e)*Xmd];
Zref(3,:) = [0, 0, rs+(p/omega_e)*Xls, 0, 0, 0, 0];
Zref(4,:) = [(p/omega_e)*Xmq, 0, 0, rkq1+(p/omega_e)*Xkq1, (p/omega_e)*Xmq, 0, 0];
Zref(5,:) = [(p/omega_e)*Xmq, 0, 0, (p/omega_e)*Xmq, rkq2+(p/omega_e)*Xkq2, 0, 0];
Zref(6,:) = [0, (Xmd/rfd)*(p/omega_e)*Xmd, 0, 0, 0, (Xmd/rfd)*(rfd+(p/omega_e)*Xfd), (Xmd/rfd)*(p/omega_e)*Xmd];
Zref(7,:) = [0, (p/omega_e)*Xmd, 0, 0, 0, (p/omega_e)*Xmd, rkd+(p/omega_e)*Xkd];

Xref(1,:) = [Xq, 0, 0, Xmq, Xmq, 0, 0];
Xref(2,:) = [0, Xd, 0, 0, 0, Xmd, Xmd];
Xref(3,:) = [0, 0, Xls, 0, 0, 0, 0];
Xref(4,:) = [Xmq, 0, 0, Xkq1, Xmq, 0, 0];
Xref(5,:) = [Xmq, 0, 0, Xmq, Xkq2, 0, 0];
Xref(6,:) = [0, Xmd, 0, 0, 0, Xfd, Xmd];
Xref(7,:) = [0, Xmd, 0, 0, 0, Xmd, Xkd];

Iref = inv(Zref)*Vref;

PhiRef = Xref*Iref;

Iqdos(1,1) = Iref(1,1);
Iqdos(2,1) = Iref(2,1);
Iqdos(3,1) = Iref(3,1);

Iabcs = inv(Ks)*Iqdos;

%% Now to enter our loop and fill in the vectors we want to plot

while t <= tmax
    
    %% We want vqs, vds, iqs, ids, P and the mechanical speeed as a function of time
    
    time(index) = t;
    Vas(index) = Vabcs(1,1);
    Vbs(index) = Vabcs(2,1);
    Vcs(index) = Vabcs(3,1);
    
    Vqs(index) = Vqdos(1,1);
    Vds(index) = Vqdos(2,1);
    Vos(index) = Vqdos(3,1);
    
    Iqs(index) = Iqdos(1,1);
    Ids(index) = Iqdos(2,1);
    Ios(index) = Iqdos(3,1);   
    
    Ias(index) = Iabcs(1,1);
    Ibs(index) = Iabcs(2,1);
    Ics(index) = Iabcs(3,1);
    
    P(index) = (Vabcs(1,1)*Iabcs(1,1)+Vabcs(2,1)*Iabcs(2,1)+Vabcs(3,1)*Iabcs(3,1));
    
    Te(index) = abs((3/2)*(p/2)*(1/omega_e)*(PhiRef(2,1)*Iref(1,1)-PhiRef(1,1)*Iref(2,1)));
    
    %% Now to take our next step through our iteration
    
    t = t + h;
    theta = theta+omega_r*h;
    index = index+1;
    
    %% And recalculate Vabcs, Iabcs, Vqdos, and Iqdos, as well as Ks
    
    Ks(1,:) = [cos(theta), cos(theta-2*pi/3), cos(theta+2*pi/3)];
    Ks(2,:) = [sin(theta), sin(theta-2*pi/3), sin(theta+2*pi/3)];   
    Ks(3,:) = [1/2, 1/2, 1/2];
    Ks = (2/3)*Ks;              %% Reference frame Matrix


    %% Now to define initial values for Vabc, Iabc, Vqdo, and Iqdo. 

    Vabcs(1,1) = sqrt(2)*Vs*cos(omega_e*t);
    Vabcs(2,1) = sqrt(2)*Vs*cos(omega_e*t-2*pi/3);
    Vabcs(3,1) = sqrt(2)*Vs*cos(omega_e*t+2*pi/3);

    Vqdos = Ks * Vabcs;
    Vkq1 = 0;
    Vkq2 = 0;
    exfd = 10;
    Vkd = 0;

    Vref(1,1) = Vqdos(1,1);        
    Vref(2,1) = Vqdos(2,1);
    Vref(3,1) = Vqdos(3,1);
    Vref(4,1) = Vkq1;
    Vref(5,1) = Vkq2;
    Vref(6,1) = exfd;
    Vref(7,1) = Vkd;

    %% Need to build this matrix to find the currents

    Zref(1,:) = [rs+(p/omega_e)*Xq, (omega_r/omega_e)*Xd, 0, (p/omega_e)*Xmq, (p/omega_e)*Xmq, (omega_r/omega_e)*Xmd, (omega_r/omega_e)*Xmd];
    Zref(2,:) = [-(omega_r/omega_e)*Xq, rs+(p/omega_e)*Xd, 0, -(omega_r/omega_e)*Xmq, -(omega_r/omega_e)*Xmq, (p/omega_e)*Xmd, (p/omega_e)*Xmd];
    Zref(3,:) = [0, 0, rs+(p/omega_e)*Xls, 0, 0, 0, 0];
    Zref(4,:) = [(p/omega_e)*Xmq, 0, 0, rkq1+(p/omega_e)*Xkq1, (p/omega_e)*Xmq, 0, 0];
    Zref(5,:) = [(p/omega_e)*Xmq, 0, 0, (p/omega_e)*Xmq, rkq2+(p/omega_e)*Xkq2, 0, 0];
    Zref(6,:) = [0, (Xmd/rfd)*(p/omega_e)*Xmd, 0, 0, 0, (Xmd/rfd)*(rfd+(p/omega_e)*Xfd), (Xmd/rfd)*(p/omega_e)*Xmd];
    Zref(7,:) = [0, (p/omega_e)*Xmd, 0, 0, 0, (p/omega_e)*Xmd, rkd+(p/omega_e)*Xkd];

    Xref(1,:) = [Xq, 0, 0, Xmq, Xmq, 0, 0];
    Xref(2,:) = [0, Xd, 0, 0, 0, Xmd, Xmd];
    Xref(3,:) = [0, 0, Xls, 0, 0, 0, 0];
    Xref(4,:) = [Xmq, 0, 0, Xkq1, Xmq, 0, 0];
    Xref(5,:) = [Xmq, 0, 0, Xmq, Xkq2, 0, 0];
    Xref(6,:) = [0, Xmd, 0, 0, 0, Xfd, Xmd];
    Xref(7,:) = [0, Xmd, 0, 0, 0, Xmd, Xkd];

    Iref = inv(Zref)*Vref;

    PhiRef = Xref*Iref;

    Iqdos(1,1) = Iref(1,1);
    Iqdos(2,1) = Iref(2,1);
    Iqdos(3,1) = Iref(3,1);

    Iabcs = inv(Ks)*Iqdos;
end


figure(1),plot(time,Ias,time,Ibs,time,Ics)

figure(2),plot(time,Iqs,time,Ids,time,Ios)

figure(3),plot(time,Vqs,time,Vds,time,Vos)

figure(4),plot(time,Vas,time,Vbs,time,Vcs)

figure(5),plot(time,P)
axis([0 0.1 0 60])

figure(6),plot(time,Te)

figure(7), plot(time,Vas,time,Ias)