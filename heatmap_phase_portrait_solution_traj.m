%%%code in this script is produced with comments explaining what the line of code next to or below it does
%%%code is explained for phase portrait array solution trajectories figures
%%%code for 1,2,3 or 4 stable steady states needs to be selected and the remaining commented out 

%pre-setting the font, figure size and fontsizes
fs=9;fn='Helvetica';wd=8;ht=7;
%setting time range to integrate over
tspan=linspace(0,100,5000);
%parameter values for ODEs that are being fixed
DegradationStrength=1;n=4;theta_a1=0.5;theta_a2=0.5;theta_b1=0.5;theta_b2=0.5;
%output directory for figures 
folder = 'U:\PhD\UoB\Figures\Heatmaps';
%ode45 tolerances
ode_options=odeset('RelTol',1e-10,'AbsTol',1e-12); %changing tolerances ODE45

%creating figure
fig1=figure;
%setting matrix row value to zero
matrix_row=0;

%setting parameter sets to get 1-4 stable steady states (sss)
% % ActivationStrength=1;InhibitionStrength=1;lambda=0.25;%1sss
% % ActivationStrength=0.5;InhibitionStrength=1;lambda=1;%2sss
% % ActivationStrength=1;InhibitionStrength=1;lambda=1;%3sss
ActivationStrength=1;InhibitionStrength=0;lambda=1;%4sss
%pre-setting matrix size to speed up computations
M1=zeros(150,5);
%ODEs
f = @(t,x) [lambda*ActivationStrength*x(1)^n./(theta_a1^n+x(1)^n)+lambda*InhibitionStrength*theta_b1^n./(theta_b1^n+x(2)^n)-DegradationStrength*x(1);...
            lambda*ActivationStrength*x(2)^n./(theta_a2^n+x(2)^n)+lambda*InhibitionStrength*theta_b2^n./(theta_b2^n+x(1)^n)-DegradationStrength*x(2)];
%symbolic variables for protein levels
syms x1 x2;
%ODEs
f_sym = [lambda*ActivationStrength*x1^n./(theta_a1^n+x1^n)+lambda*InhibitionStrength*theta_b1^n./(theta_b1^n+x2^n)-DegradationStrength*x1;lambda*ActivationStrength*x2^n./(theta_a2^n+x2^n)+lambda*InhibitionStrength*theta_b2^n./(theta_b2^n+x1^n)-DegradationStrength*x2];
%variables for jacobian matrix
v_sym=[x1,x2];
%calculating jacbian with respect to variables x1 & x2
jac=jacobian(f_sym,v_sym);
    
	for i=0:0.4:4%initial conditions on x axis
		%display where the computation is up to in command window - good for long computations to see where you are up to
        fprintf('Running first subplot with i=%.2f.\n',i);
        for j=0:0.4:4%initial conditions on x axis
			%increasing the value by 1 in each new loop - moves the row
			%values up by one each loop
            matrix_row=matrix_row+1;
			%using ode45
            [t,x_num]=ode45(f,tspan,[i,j],ode_options); %solving ODEs with ics
			%calculated steady state values
            x1_ss=x_num(5000,1);x2_ss=x_num(5000,2);
			%rounded calculated steady state values
            x1_ss_b=round(x1_ss,3);x2_ss_b=round(x2_ss,3); %rounding steady state(ss) position to see unique ss
            
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
            
            %matrix of ics, ss positions and the stability
            M1(matrix_row,:) = [i j x1_ss_b x2_ss_b stability];
            
			%axes limits & box around figure
            grid on;hold on;box on;
			
			%plotting solution trajectory if stability is stable
            if stability == 1
% % 				 plot(x_num(:,1),x_num(:,2),'Color',[0.1,0.9,1]);%1ss
% %                  plot(x_num(:,1),x_num(:,2),'Color',[0,0.5,1]);%2ss
% %                  plot(x_num(:,1),x_num(:,2),'Color',[0.4,0.1,1]);%3ss
                 plot(x_num(:,1),x_num(:,2),'m');%4ss
			else
				point = 0;
            end
        end
	end

%selecting a column in M1 matrix
col_stable=M1(:,5);
%new sub-matrix M1_b is a submatrix of M1 with stability value = 1 (stable stedy states)
M1_b=M1(col_stable==1,:);
%extracting unique stable steady states 
M2 = unique(M1_b(:,[3 4]),'rows');

%plotting stable steady states
plot(M2(:,1),M2(:,2),'o','MarkerEdgeColor','black', 'MarkerFaceColor',[0.49 1 0.63],...
    'MarkerSize',4);
%removing axis labels
set(gca,'YTickLabel',[],'XTickLabel',[]);
%setting figure size
set(gcf,'Units','centimeters','Position',[0 0 wd ht],'PaperUnits','centimeters','PaperSize',[wd ht]);
%setting axis limits
xlim([0 4]);ylim([0 4]);hold off;

%saving produced figure to output directory with specified name and file extenstion
% % epsFileName1 = '1_steady_state.eps';fullFileName =fullfile(folder, epsFileName1);print(fig1,fullFileName,'-depsc');%1sss
% % epsFileName1 = '2_steady_states.eps';fullFileName =fullfile(folder, epsFileName1);print(fig1,fullFileName,'-depsc');%2sss
% % epsFileName1 = '3_steady_states.eps';fullFileName =fullfile(folder, epsFileName1);print(fig1,fullFileName,'-depsc');%3sss
% % epsFileName1 = '4_steady_states.eps';fullFileName =fullfile(folder, epsFileName1);print(fig1,fullFileName,'-depsc');%4sss
