function stabilitydiagram(tminutes);
%stabilitydiagram(240); %integrate trajectories for at most 240 minutes
thours = tminutes/60;
figure(3); clf; hold on;
options = odeset('RelTol',10^-6,'AbsTol',10^-9);
mlist = linspace(200,500,21);
Tlist = linspace(350,560,21);
for(i = 1:length(mlist))
    for(j=1:length(Tlist))
        [t,x] = ode15s(@odefun,[0,thours],[mlist(i),Tlist(j)],options);
        m = x(:,1);
        T = x(:,2);
        plot(m(end),T(end),'or'); %red stop
        plot(m(1),T(1),'og'); %green go
        plot(m,T,'-k'); %black trajectory
    end
end
title('stability diagram at flow shutoff')
xlabel('NH_4NO_3 holdup mass, kg')
ylabel('reactor temperature, degF')
xlim([0,500])
ylim([0,600])
grid on;
hold off;
return
function derivs = odefun(t,x)
m = x(1); %kg NH4NO3
T = x(2); %degF
k = 2.28e19 *exp(-44367/(T+459.67));
UA = 0.01073; %MJ/(h degF)
Ta = 100; %degF
cpa = 0.8/1000; %MJ/(kg NH4NO3 degF)
dHrxn = -0.740; %MJ/kg NH4NO3
dmdt = -k*m;
%dTdt = (-k*dHrxn/cpa) + (UA*(Ta-T)/(m*cpa)) + (T*k); %slightly wrong
dTdt = (-k*dHrxn/cpa) + (UA*(Ta-T)/(m*cpa)); %correct
derivs = [dmdt;dTdt];
return;