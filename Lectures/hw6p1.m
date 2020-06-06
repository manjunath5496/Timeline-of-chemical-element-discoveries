function [m,T]=hw6p1(tminutes);
%[m,T]=hw6p1(5); %run simulation for 5 minutes, crashes at ~4.5 minutes
thours = tminutes/60;
x0 = [378.5;510];
options = odeset('RelTol',10^-6,'AbsTol',10^-9);
[t,x] = ode15s(@odefun,[0,thours],x0,options);
m = x(:,1);
T = x(:,2);
figure(1);
plot(t*60,m);
title('NH_4NO_3 holdup');
xlabel('time/minutes since feed shutoff')
ylabel('m')
figure(2);
plot(t*60,T);
title('Temperature in Reactor');
xlabel('time/minutes since feed shutoff')
ylabel('T/degF')
return;
function derivs = odefun(t,x)
m = x(1); %kg NH4NO3
T = x(2); %degF
k = 2.28e19 *exp(-44367*(0.99)/(T+459.67));
UA = 0.01073; %MJ/(h degF)
Ta = 100; %degF
cpa = 0.8/1000; %MJ/(kg NH4NO3 degF)
dHrxn = -0.740; %MJ/kg NH4NO3
dmdt = -k*m;
%dTdt = (-k*dHrxn/cpa) + (UA*(Ta-T)/(m*cpa)) + (T*k); %slightly wrong
dTdt = (-k*dHrxn/cpa) + (UA*(Ta-T)/(m*cpa)); %correct
derivs = [dmdt;dTdt];
return;

