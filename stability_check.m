%%%code in this script is produced with comments explaining what the line of code next to or below it does

%fixed parameter values for ODEs
DegradationStrength=1;n=4;theta_a1=0.5;theta_a2=0.5;theta_b1=0.5;theta_b2=0.5;
%variable parameters in ODEs
ActivationStrength=1;InhibitionStrength=0.25;

%lambda, if known
lambda=0.961;
%ATP value and then lambda calculated from this
% % ATP=1500;%ATP values
% % s1=0.5;s2=0.0033;s3=-5;s4=0.5;%paramter values within lambda function
% % l= @(ATP) s1*tanh(s2*ATP+s3)+s4;%lambda

%calculated steady state
steadystate=[.3609751191 2.434466763];
%symbolic variables for protein levels
syms x1 x2;
%ODEs
ODEs = [lambda*ActivationStrength*x1^n./(theta_a1^n+x1^n)+lambda*InhibitionStrength*theta_b1^n./(theta_b1^n+x2^n)-DegradationStrength*x1;lambda*ActivationStrength*x2^n./(theta_a2^n+x2^n)+lambda*InhibitionStrength*theta_b2^n./(theta_b2^n+x1^n)-DegradationStrength*x2];
%alternative ODEs when calculating using ATP
% % ODEs=[l(ATP)*ActivationStrength*x1^n./(theta_a1^n+x1^n)+l(ATP)*InhibitionStrength*theta_b1^n./(theta_b1^n+x2^n)-DegradationStrength*x1;l(ATP)*ActivationStrength*x2^n./(theta_a2^n+x2^n)+l(ATP)*InhibitionStrength*theta_b2^n./(theta_b2^n+x1^n)-DegradationStrength*x2];

%variables for jacobian
variables=[x1,x2];
%calculating jacbian with respect to variables x1 & x2
jac=jacobian(ODEs,variables); 
%substitute steady state into Jacobian matrix
sub=subs(jac, [x1 x2], steadystate);
%calculate eigenvalues
eigenvalues = eig(sub);
%sign of both eigenvalues in command window
eigen_sign=([sign(eigenvalues(1)),sign(eigenvalues(2))]);
%display sign of both eigenvalues in command window
disp(eigen_sign);
