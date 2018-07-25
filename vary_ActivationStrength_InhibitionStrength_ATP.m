%declaring a function named vary_ActivationStrength_InhibitionStrength_ATP that accepts inputs X,ActivationStrength,InhibitionStrength & ATP
% and returns outputs in F
function F = vary_ActivationStrength_InhibitionStrength_ATP(X,ActivationStrength,InhibitionStrength,ATP)

%setting variable names for protein levels in ODEs
x1=X(1);
x2=X(2);

%pre-setting fixed parameter values
DegradationStrength=1;
n=4;theta_a1=0.5;theta_a2=0.5;theta_b1=0.5;theta_b2=0.5;

%setting lambda function
s1=0.5;s2=0.0033;s3=-5;s4=0.5;
l= @(ATP) s1*tanh(s2*ATP+s3)+s4;

%ODEs for x1 and x2 respectively
x1dot=l(ATP)*ActivationStrength*x1^n./(theta_a1^n+x1^n)+l(ATP)*InhibitionStrength*theta_b1^n./(theta_b1^n+x2^n)-DegradationStrength*x1;
x2dot=l(ATP)*ActivationStrength*x2^n./(theta_a2^n+x2^n)+l(ATP)*InhibitionStrength*theta_b2^n./(theta_b2^n+x1^n)-DegradationStrength*x2;

%function of both ODEs together to make an ODE system of equations
F=[x1dot;x2dot];