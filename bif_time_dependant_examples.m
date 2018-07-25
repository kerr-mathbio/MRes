%%%code in this script is produced with comments explaining what the line of code next to or below it does
%%%bifurcation diagram time-dependant figures

%pre-setting the font, figure size and fontsizes
fn='Helvetica';wd=8;ht=7;fs_labels=12;fs_axis=10;
%parameter values for ODEs that are being fixed
ActivationStrength=1;InhibitionStrength=1;DegradationStrength=1;n=4;theta_a1=0.5;theta_a2=0.5;theta_b1=0.5;theta_b2=0.5;
%parameter values for sigmoidal curve
s1=0.5;s2=0.0033;s3=-5;s4=0.5;
%ode45 tolerances
options = odeset('RelTol',1e-10,'AbsTol',1e-12);
%time range to integrate over
tspan=linspace(0,100,2000);
%output directory for figures
folder= 'U:\PhD\UoB\Figures';

%ATP levels to plot diagrams for 
ATP=2200;%1200;1520;
%lambda
l= @(ATP) s1*tanh(s2*ATP+s3)+s4;%lambda([ATP])
%ODEs
odes = @(t,x) [l(ATP)*ActivationStrength*x(1)^n./(theta_a1^n+x(1)^n)+l(ATP)*InhibitionStrength*theta_b1^n./(theta_b1^n+x(2)^n)-DegradationStrength*x(1);...
                l(ATP)*ActivationStrength*x(2)^n./(theta_a2^n+x(2)^n)+l(ATP)*InhibitionStrength*theta_b2^n./(theta_b2^n+x(1)^n)-DegradationStrength*x(2)];
%symbolic variables for protein levels
syms x1 x2;
%ODEs
f_sym = [l(ATP)*ActivationStrength*x1^n./(theta_a1^n+x1^n)+l(ATP)*InhibitionStrength*theta_b1^n./(theta_b1^n+x2^n)-DegradationStrength*x1;l(ATP)*ActivationStrength*x2^n./(theta_a2^n+x2^n)+l(ATP)*InhibitionStrength*theta_b2^n./(theta_b2^n+x1^n)-DegradationStrength*x2];
% variables for jacobian matrix
v_sym=[x1,x2];
%calculating jacobian with respect to variables x1 & x2
jac=jacobian(f_sym,v_sym);

%creating figure
fig1=figure(1);clf;hold on;
    for i=0:0.5:5%initial conditions on x axis
        for j=0:0.5:5%initial conditions on y axis
            %initial conditions pairing
			ics=[i,j];
			%using ode45
            [t,x_num]=ode45(odes,tspan,ics,options);
			%calculated steady state values
			x1_ss=x_num(2000,1);x2_ss=x_num(2000,2);
			%rounded calculated steady state values
			x1_ss2=round(x1_ss,3);x2_ss2=round(x2_ss,3);

			%subs. in steady state values to jacobian
			sub=subs(jac, [x1 x2], [x1_ss x2_ss]); %subs. in ss values from original ics
			%calc eigenvlaues
			eigen = eig(sub); %calc eigenvlaues of matrix 'sub'
			%calculate the sign of each eigenvalue
			eigenvalue_1=sign(eigen(1));eigenvalue_2=sign(eigen(2));
			
			%testing if the steady state is stable or unstable
			if (eigenvalue_1 < 0) && (eigenvalue_2 < 0)
				stability = 1;
			else 
				stability = -1;
			end

			%x and y axis limits + box and grid for figure
            xlim([0 15]);ylim([0 5]);grid on;box on;
			%plotting time-dependant solution if end steady state is stable
			if stability == 1
				 plot(t,x_num(:,2),'b');
			else
				point = 0;
			end
            ax = gca;ax.FontSize=fs_axis;ax.FontName=fn;ax.TickDir = 'out';%changing x and y axes properties
            fig = gcf;fig.Units='centimeters';fig.Position=[0 0 wd ht];fig.PaperUnits='centimeters';fig.PaperSize=[wd ht];%setting figure size
        end
    end
hold off;
%saving produced figure to output directory with specified name and file extenstion
pngFileName = sprintf('bif_ex_ATP=%d.eps',ATP);fullFileName=fullfile(folder,pngFileName);print(fig1,fullFileName,'-depsc');