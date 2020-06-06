function [t,n,X,V,xi1,xi2,xi3]= hwk1prob1(k1,k2,k3,tmax,reltol,abstol);
%[t,n,X,V,xi1,xi2,xi3] = hwk1prob1(k1,k2,k3,tmax,reltol,abstol);
%EXAMPLE
%[t,n,X,V,xi1,xi2,xi3] = hwk1prob1(10^3,10^4,10^5,(10)/10^3,1e-6,1e-9);
%Calculate airbag dynamics with an initial charge of NaN3,NaNO3,and SiO2
%
%Plot extents of reaction, airbag volume, species amounts, and conversion
%vs time
%
%k1: 10^3 mol/s
%k2: 10^4 L/(mol s)
%k3: 10^5 L/s
%tmax: 10/10^3 s (this is 10/k1, 10 milliseconds)
%tmax: 0.1/10^3 s (this is 0.1 milliseconds)

%PREPROCESSING-------------------------------------------------------------
%initial amounts of each species stored in struct "n0"
n0.NaN3  = 150/65;
n0.NaNO3 =(n0.NaN3)/5;
n0.SiO2  = 6*(n0.NaN3);
n0.Na    = 0;
n0.N2    = 0;
n0.Na2O  = 0;
n0.glass = 0;

%problem parameters (constants) stored in struct "param"
%so they are available in the ODE solver
param.V0     = 70/1000; %L
param.VN     = 22;      %L/mol
param.csound = 300 * 100 / 10; %L^(1/3)/s; 10 cm/s = ( (1L)^(1/3) )/s
param.n0     = n0;      %mol
param.k1     = k1;      %mol/s
param.k2     = k2;      %k2: L/(mol s)
param.k3     = k3;      %k3: L/s
param.reltol = reltol;
param.abstol = abstol;

%set initial condition
%x0=     [xi1_0,xi2_0,xi3_0,      V0]';
x0 =     [    0,    0,    0,param.V0]';

%INTEGRATION---------------------------------------------------------------

%call ode15s integrator with tolerance options;
options = odeset('RelTol',reltol,'AbsTol',abstol);
[t,x] = ode15s(@reactionfun,[0,tmax],x0,options,param);

%POSTPROCESSING------------------------------------------------------------
%rename dependent variables
xi1 = x(:,1);
xi2 = x(:,2);
xi3 = x(:,3);
V   = x(:,4);

%calculate moles of each species at all times
n.NaN3   = n0.NaN3  -  2 * xi1;
n.NaNO3  = n0.NaNO3 -  2 * xi2;
n.SiO2   = n0.SiO2  - 10 * xi3;
n.Na     = 2 * xi1  - 10 * xi2;
n.N2    =  3 * xi1  +      xi2;
n.Na2O  =  6 * xi2  -      xi3;
n.glass =      xi3;

%calculate conversion at all times
%1 - (final NaN3 / initial NaN3)
X.NaN3 =  1-(n.NaN3)/(n0.NaN3);
%(N atoms in Nitrogen)/(initial N atoms)
X.N2   =  (2*n.N2)/(3*n0.NaN3 + n0.NaNO3);

%plot extents of reaction vs time
figure(1);
plot(t,xi1,t,xi2,t,xi3);
xlabel('time/s');
ylabel('moles');
legend('xi1','xi2','xi3')

%plot airbag volume vs time
figure(2);
plot(t,V);
xlabel('time/s')
ylabel('V/L')

%plot species amounts vs time
figure(3);
plot(t,n.NaN3,...
    t,n.NaNO3,...
    t,n.SiO2,...
    t,n.Na,...
    t,n.N2,...
    t,n.Na2O,...
    t,n.glass);
legend('NaN3','NaNO3','SiO2','Na','N2','Na2O','glass');
xlabel('time/s')
ylabel('moles')

%plot extents of reaction vs time
figure(4);
plot(t,X.NaN3,...
    t,X.N2);
xlabel('time/s')
ylabel('conversion')
legend('NaN3','N2')
return;
%--------------------------------------------------------------------------
function derivs = reactionfun(t,x,param)
%retrieve constants from param
k1 = param.k1;
k2 = param.k2;
k3 = param.k3;
V0 = param.V0;
VN = param.VN;
csound = param.csound;
n0_NaN3 = param.n0.NaN3;
n0_NaNO3 = param.n0.NaNO3;
n0_SiO2 = param.n0.SiO2;

%rename dependent variables for clarity
xi1 = x(1);
xi2 = x(2);
xi3 = x(3);
V   = x(4);

%calculate moles of each species in terms of extents of reaction
n.NaN3   = n0_NaN3  -  2 * xi1;
n.NaNO3  = n0_NaNO3 -  2 * xi2;
n.SiO2   = n0_SiO2  - 10 * xi3;
n.Na     = 2 * xi1  - 10 * xi2;
n.N2    =  3 * xi1  +      xi2;
n.Na2O  =  6 * xi2  -      xi3;
n.glass =      xi3;

%calculate rate of change of dependent variables in terms of species
%concentrations and retrieved parameters
%if(n.NaN3 > small fraction of initial amount present) react normally
if(n.NaN3 > 0.001*n0_NaN3)
    dxi1dt = k1;
else %shut down zero-order dynamics
    dxi1dt = 0;
end
dxi2dt = k2 * (n.Na/V) * (n.NaNO3/V) * V;
dxi3dt = k3 * (n.Na2O/V);

if(V < V0 + VN*n.N2)
    dVdt = (36*pi*V^2)^(1/3)*csound;
else
    dVdt = VN * (3*dxi1dt + dxi2dt);
end

%gather derivatives in one column vector for output
derivs = [dxi1dt; dxi2dt; dxi3dt; dVdt];

return;
%--------------------------------------------------------------------------