function [k,k0]= hwk1prob2(C2H40_1,C2H40_2,C2H40_3);
%hwk1prob2(C2H40_1,C2H40_2,C2H40_3)
%hwk1prob2(6.7e-4,4e-4,1.33e-4)
%
%Calculate k for reaction of vinyl radical (C2H3) with ethene (C2H4) using
%three sets of experimental data (time vs signal).
%Compare model predictions with experimental data.
%
%Initial C2H4 concentrations for each experimental data set
%[C2H4]0,1: 6.7x10^-4 M
%[C2H4]0,2: 4.0x10^-4 M
%[C2H4]0,3: 1.33x10^-4 M
%
%k has units of L/mol-s
%k0 has units of 1/s

%initial [C2H4] for each set of experimental data
C2H40=[C2H40_1; C2H40_2; C2H40_3];
    
%FIT EXPERIMENTAL DATA (for Bn, An and taun)
%Load experimental data from text files
%The first column of each text file is time, second column is signal.
    for i=1:3
        dataset = ['vinylethene' int2str(i) '.txt'];
        data=load(dataset);
        time=data(:,1);
        Sn=data(:,2);

%define model specify "time" as the independent variable
%set initial guesses for Bn, An and taun

        initial_guess=[.0003,.0003,.00012];
        model=fittype('Bn+An*exp(-time/taun)','ind','time');
        options=fitoptions(model);
        set(options,'StartPoint',initial_guess);
        [fresult] = fit(time,Sn,model,options);
        An(i)=fresult.An;
        Bn(i)=fresult.Bn;
        taun(i)=fresult.taun;
        inv_taun(i,1)=1/taun(i);
    end

%Plot the experimental 1/taun vs [C2H4]0 (three data points).
figure(1);
plot(C2H40,inv_taun,'x');
xlabel('[C2H4]0 (M)')
ylabel('1/taun (1/s)')
hold on

%Fit the 1/taun vs [C2H4]0 data to a straight line to determine k and r0.
%Provide intial r0 and k values.
%Define the model and set C2H40 as the independent variable.
initial_guess2=[100, 10^7];
model2=fittype('k0+k*C2H40','ind','C2H40');
options=fitoptions(model2);
set(options,'StartPoint',initial_guess2);
[fresult]=fit(C2H40,inv_taun,model2,options);
k=fresult.k;  %k has units of L/mol-s
k0=fresult.k0;  %k0 has units of 1/s

%Calculate 1/taun for a range of [C2H4]0 values using the k0 and k
%parameters found above.
for j=1:10
    C2H40_calc(j)=C2H40_1*j/10;
    inv_tau_calc(j)=k0+k*C2H40_calc(j);
end

%Add the calculated 1/taun data as a line on the experimental 1/taun vs
%[C2H4]0 plot.  
plot(C2H40_calc,inv_tau_calc);
legend('exp','fit')
hold off

%Compare experimental data with calculated signal.
load vinylethene1.txt
time_1=vinylethene1(:,1);
Sn_1=vinylethene1(:,2);

load vinylethene2.txt
time_2=vinylethene2(:,1);
Sn_2=vinylethene2(:,2);

load vinylethene3.txt
time_3=vinylethene3(:,1);
Sn_3=vinylethene3(:,2);

%Calculate tau_n for each initial [C2H4]_n using k0 and k.
for i=1:3
    tau(i)=1/(k0+k*C2H40(i));
end

%Calculate signal from Bn, An and taun parameters.
signal_calc1=Bn(1)+An(1)*exp(-time_1/tau(1));
signal_calc2=Bn(2)+An(2)*exp(-time_2/tau(2));
signal_calc3=Bn(3)+An(3)*exp(-time_3/tau(3));

%Plot the experimental and calculated signal vs time for each set of data.
figure(2);
plot(time_1,Sn_1,'r.');
xlabel('time(s)');
ylabel('signal');
title('Vinylethene 1')
hold on
plot(time_1,signal_calc1);
legend('exp','calc')
hold off

figure (3);
plot(time_2,Sn_2,'r.');
xlabel('time(s)');
ylabel('signal');
title('Vinylethene 2')
hold on
plot(time_2,signal_calc2);
legend('exp','calc')
hold off

figure(4);
plot(time_3,Sn_3,'r.');
xlabel('time(s)');
ylabel('signal');
title('Vinylethene 3')
hold on
plot(time_3,signal_calc3);
legend('exp','calc')
hold off

return;