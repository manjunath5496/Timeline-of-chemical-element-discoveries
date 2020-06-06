function [FB1rctr,FB3rctr] = hw3prob3();
%[FB1rctr,FB3rctr] = hw3prob3
param.k1  = 1.82e-4*3600; %1/h
param.km1 = 4.49e-4*3600; %1/M 1/h

%WARNING: the values stored in q, ca, and cb variables change as the
%program executes

% part A
param.V = 15; %L
q0 = 0.25; %L/h
ca0 = 0.25; %M
cb0 = 0;
[ca,cb] = reactorX(ca0,cb0,q0,param);
FB1rctr = q0 * cb; %mol/h

% part B
FB3rctr = 0; %accumulator for B production, mol/h
param.V = 5; %L
q0 = 0.25; %L/h
ca0 = 0.25; %M
cb0 = 0;

for(i = 1:2) %reactor1, separator1, reactor2, separator2
    [ca,cb] = reactorX(ca0,cb0,q0,param);
    %rename reactor effluent as separator input
    ca0 = ca;
    cb0 = cb;
    q0 = q0;

    [ca1,cb1,q1,FBE] = separatorX(ca0,cb0,q0);

    FB3rctr = FB3rctr + FBE;

    %rename separator effluent as reactor input
    ca0 = ca1;
    cb0 = cb1;
    q0 = q1;
end

%reactor3
[ca,cb] = reactorX(ca0,cb0,q0,param);

FB3rctr = FB3rctr + cb*q0;

return;

function [ca,cb] = reactorX(ca0,cb0,q0,param)
    k1V =  param.k1  * param.V; %L/h
    km1V = param.km1 * param.V; %L^2/(mol h)
    
    %quadratic equation in cb: coefficients of polynomial
    a = km1V;
    b = q0 + k1V;
    c = -q0*cb0 - k1V*ca0 - k1V*cb0;
    
    cb = (-b + sqrt(b^2 - 4*a*c))/(2*a);
    ca = ca0 + cb0 - cb;

return;

function [ca1,cb1,q1,FBE]=separatorX(ca0,cb0,q0)
    q1overq2 = (ca0 + 0.5*cb0)/(1.5*cb0);
    q2 = q0/(1 + q1overq2);
    q1 = q0 - q2;

    FBE = 0.75 * cb0 * q0; %B recovered from bottom stream
    
    ca1 = q0/q1 * ca0;
    cb1 = 0.25 * q0/q1 * cb0;

return;