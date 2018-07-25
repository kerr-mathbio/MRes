%%%code in this script is produced with comments explaining what the line of code next to or below it does
%%%code is explained for bifurcation b diagrams and then are very similar for a and k bifurcation diagrams

%pre-setting the font, figure size and fontsizes
fn='Helvetica';wd=8;ht=7;fs_labels=10;fs_axis=9;
%parameter values for ODEs that are being fixed
n=4;theta_a1=0.5;theta_a2=0.5;theta_b1=0.5;theta_b2=0.5;
%parameter values for sigmoidal curve
s1=0.5;s2=0.0033;s3=-5;s4=0.5;
%fsolve tolerances
options = optimoptions('fsolve','FunctionTolerance',1e-11,'OptimalityTolerance',1e-11,'Display','off'); 
%output directory for figure
folder= 'U:\PhD\UoB\Figures\Bifurcation_Diagrams\Combined_Param';


%%bifurcation diagram with b
%setting fixed parameter values for these bif diagram
ActivationStrength=1;DegradationStrength=1;
%pre-setting matrix size to speed up computations and setting inital row in matrix as zero
matrix_row1 = 0;ss_1a=zeros(400000,10);
for InhibitionStrength=[0.25,0.5,0.75,1.5]%b values to scan through
%creating bifurcation figure
fig1=figure('Name','Bifurcation_b');
%starting time of loop
clock_start = datestr(now,'HH:MM:SS');
	for ATP = 320:10:2760%ATP values to scan through and step size
		%display where the computation is up to in command window - good for long computations to see where you are up to
		fprintf('Running BD b with b=%.2f & [ATP]=%d at %s.\n',InhibitionStrength,ATP,datestr(now,'HH:MM:SS'));
		%lambda
		l1= @(ATP) s1*tanh(s2*ATP+s3)+s4;%lambda([ATP])
		%function to be used by fsolve 
		fhandle=@(X)bif_b_y(X,InhibitionStrength,ATP);
		
		%symbolic variables for protein levels
		syms x1 x2;
		%ODEs
		f1 = [l1(ATP)*ActivationStrength*x1^n./(theta_a1^n+x1^n)+l1(ATP)*InhibitionStrength*theta_b1^n./(theta_b1^n+x2^n)-DegradationStrength*x1;...
			l1(ATP)*ActivationStrength*x2^n./(theta_a2^n+x2^n)+l1(ATP)*InhibitionStrength*theta_b2^n./(theta_b2^n+x1^n)-DegradationStrength*x2];
		% variables for jacobian matrix
		v1=[x1,x2];
		%calculating jacbian with respect to variables x1 & x2
		jac1=jacobian(f1,v1); 
		
		for i=0:0.25:5%initial conditions on x axis
			for j=0:0.25:5%initial conditions on y axis
			%increasing matrix row by 1 in each loop
			matrix_row1=matrix_row1+1;
			
			%initial conditions pairing
			X0 = [i,j];
			%using fsolve
			X = fsolve(fhandle,X0,options);
			%calculated steady state values
			x1_ss1=X(1);x2_ss1=X(2);
			%rounded calculated steady state values
			x1_ss_1a=round(X(1),3);x2_ss_1a=round(X(2),3);
            
			%subs. in steady state values to jacobian
			sub1=subs(jac1, [x1 x2], [x1_ss1 x2_ss1]);
			%calc eigenvlaues
			eigen1 = eig(sub1);
			%calculate the sign of each eigenvalue
			eigenvalue_1_a=sign(eigen1(1));eigenvalue_2_a=sign(eigen1(2));
            %substituting the calculated steady state back into ODEs
			sub1b=subs(f1,[x1 x2],[x1_ss1 x2_ss1]);
			
			%testing if the steady state is stable or unstable
			if (eigenvalue_1_a < 0) && (eigenvalue_2_a < 0)
                stability = 1;
            else 
                stability = -1;
			end
            
			%matrix of b, ATP, lambda, ics, ss positions, the stability and subs2
			ss_1a(matrix_row1,:) = [InhibitionStrength ATP l1(ATP) i j x1_ss_1a x2_ss_1a stability sub1b(1) sub1b(2)];
			end
		end
	%checking if the steady state when sub into ODEs is within a tolerance of plus-minus 10^-10
	col_check=ss_1a(:,9);ss_1a2=ss_1a(col_check < 1e-10,:);
	col_check2=ss_1a2(:,10);ss_1a3=ss_1a2(col_check2 < 1e-10,:);
	col_check3=ss_1a3(:,9);ss_1a4=ss_1a3(col_check3 > -1e-10,:);
	col_check4=ss_1a4(:,10);ss_1a5=ss_1a4(col_check4 > -1e-10,:);
	
	
	col_b=ss_1a5(:,1);%selecting a column in ss_1a5 matrix
	ss_1b=ss_1a5(col_b==InhibitionStrength,:);%new sub-matrix ss_1b is a submatrix of ss_1a5 with value b=b in the first column
	
	%extracting stable and unstable sub-matrices
	col_stable1=ss_1b(:,8);ss_1c=ss_1b(col_stable1==1,:);ss_1e=unique(ss_1c(:,[1 2 6 7]),'rows');
	col_unstable1=ss_1b(:,8);ss_1d=ss_1b(col_unstable1==-1,:);ss_1f=unique(ss_1d(:,[1 2 6 7]),'rows');
	
	%plotting sub-figre 
	xlim([0 3000]);hold on;ylim([0 3]);grid on;box on;%axes limits & box around figure
	plot(ss_1e(:,2),ss_1e(:,4),'b.','MarkerSize',0.7);%plotting stable steady states
	plot(ss_1f(:,2),ss_1f(:,4),'r.','MarkerSize',0.7);%plotting unstable steady states
	hx=xlabel('$[ATP] /\mu$M','interpreter','latex');hx.FontSize=fs_labels;hx.FontName=fn;%label on x-axis
	hy=ylabel('$x_{2}$ Steady State Position','interpreter','latex');hy.FontSize=fs_labels;hy.FontName=fn;%label on y-axis
	ax = gca;ax.FontSize=fs_axis;ax.FontName=fn;ax.TickDir = 'out';ax.XTick = 0:500:3000;%changing x and y axes properties
    fig = gcf;fig.Units='centimeters';fig.Position=[0 0 wd ht];fig.PaperUnits='centimeters';fig.PaperSize=[wd ht];%setting figure size
    hold off;
	end
%saving produced figure to output directory with specified name and file extenstion
pngFileName = sprintf('bif_ATP_b=%.0f.eps',InhibitionStrength*100);fullFileName=fullfile(folder,pngFileName);print(fig1,fullFileName,'-depsc');

%end time of loop
clock_end = datestr(now,'HH:MM:SS');
%displays start and end time of loop computation in command window
fprintf('b=%.2f loop start = %s & end = %s.\n',InhibitionStrength,clock_start, clock_end);
end

%clear some of the information stored by matlab 
param={'b','a','k','y','i','j',...
    'x1','x2','X','X0','stability'};clear(param{:});

%%bifurcation diagram with a
InhibitionStrength=1;DegradationStrength=1;
matrix_row2=0;ss_2a=zeros(400000,10);
for ActivationStrength=[0.25,1,1.25,1.75,3]
    fig2=figure('Name','Bifurcation_a');
	clock_start = datestr(now,'HH:MM:SS');
	for ATP = 320:10:2760
        fprintf('Running BD a with a=%.2f & [ATP]=%d at %s.\n',ActivationStrength,ATP,datestr(now,'HH:MM:SS'));
        l2= @(ATP) s1*tanh(s2*ATP+s3)+s4;
		
		fhandle=@(X)bif_a_y(X,ActivationStrength,ATP);
		
		syms x1 x2;
		f2 = [l2(ATP)*ActivationStrength*x1^n./(theta_a1^n+x1^n)+l2(ATP)*InhibitionStrength*theta_b1^n./(theta_b1^n+x2^n)-DegradationStrength*x1;...
			l2(ATP)*ActivationStrength*x2^n./(theta_a2^n+x2^n)+l2(ATP)*InhibitionStrength*theta_b2^n./(theta_b2^n+x1^n)-DegradationStrength*x2];
		v2=[x1,x2];
		jac2=jacobian(f2,v2);
				
		for i=0:0.25:5
			for j=0:0.25:5
                matrix_row2=matrix_row2+1;

				X0 = [i,j];
                X = fsolve(fhandle,X0,options);
				x1_ss2=X(1);x2_ss2=X(2);
				x1_ss_2a=round(X(1),3);x2_ss_2a=round(X(2),3);

                sub2=subs(jac2, [x1 x2], [x1_ss2 x2_ss2]);
                eigen2 = eig(sub2);
                eigenvalue_1_2a=sign(eigen2(1));eigenvalue_2_2a=sign(eigen2(2));
				
				sub2b=subs(f2,[x1 x2],[x1_ss2 x2_ss2]);
				
                if (eigenvalue_1_2a < 0) && (eigenvalue_2_2a < 0)
                    stability = 1;
                else 
                    stability = -1;
                end

                ss_2a(matrix_row2,:) = [ActivationStrength ATP l2(ATP) i j x1_ss_2a x2_ss_2a stability sub2b(1) sub2b(2)];
			end
		end
	col_check=ss_2a(:,9);ss_2a2=ss_2a(col_check < 1e-10,:);
	col_check2=ss_2a2(:,10);ss_2a3=ss_2a2(col_check2 < 1e-10,:);
	col_check3=ss_2a3(:,9);ss_2a4=ss_2a3(col_check3 > -1e-10,:);
	col_check4=ss_2a4(:,10);ss_2a5=ss_2a4(col_check4 > -1e-10,:);
	
	col_a=ss_2a5(:,1);ss_2b=ss_2a5(col_a==ActivationStrength,:);
	
	col_stable2=ss_2b(:,8);ss_2c=ss_2b(col_stable2==1,:);ss_2e=unique(ss_2c(:,[1 2 6 7]),'rows');
	col_unstable2=ss_2b(:,8);ss_2d=ss_2b(col_unstable2==-1,:);ss_2f=unique(ss_2d(:,[1 2 6 7]),'rows');
	
	xlim([0 3000]);hold on;ylim([0 4]);grid on;box on;
	plot(ss_2e(:,2),ss_2e(:,4),'b.','MarkerSize',0.7);
	plot(ss_2f(:,2),ss_2f(:,4),'r.','MarkerSize',0.7);
	hx=xlabel('$[ATP] /\mu$M','interpreter','latex');hx.FontSize=fs_labels;hx.FontName=fn;
	hy=ylabel('$x_{2}$ Steady State Position','interpreter','latex');hy.FontSize=fs_labels;hy.FontName=fn;
	ax = gca;ax.FontSize=fs_axis;ax.FontName=fn;ax.TickDir = 'out';ax.XTick = 0:500:3000;
    fig = gcf;fig.Units='centimeters';fig.Position=[0 0 wd ht];fig.PaperUnits='centimeters';fig.PaperSize=[wd ht];
    tit1=title(sprintf('a=%.2f, b=%.2f & k=%.2f',ActivationStrength,InhibitionStrength,DegradationStrength'));
    set(tit1,'fontsize',fs_axis_title);set(tit1,'fontname',fn);
    hold off;
	end
pngFileName3 = sprintf('bif_ATP_a=%.0f.eps',ActivationStrength*100);fullFileName3=fullfile(folder,pngFileName3);print(fig2,fullFileName3,'-depsc');

clock_end = datestr(now,'HH:MM:SS');
fprintf('a=%.2f loop start = %s & end = %s.\n',ActivationStrength,clock_start, clock_end);
end

param={'a','b','k','y','i','j',...
    'x1','x2','X','X0','stability'};clear(param{:});

%%bifurcation diagram with k
ActivationStrength=1;InhibitionStrength=1;
matrix_row3=0;ss_3a=zeros(400000,10);
for DegradationStrength=[0.5,1.25,1.5,3]
    fig3=figure('Name','Bifurcation_k');
	clock_start = datestr(now,'HH:MM:SS');
	for ATP = 320:10:2760
        fprintf('Running BD k with k=%.2f & [ATP]=%d at %s.\n',DegradationStrength,ATP,datestr(now,'HH:MM:SS'));
        l3= @(ATP) s1*tanh(s2*ATP+s3)+s4;
		
		fhandle=@(X)bif_k_y(X,DegradationStrength,ATP);
		
		syms x1 x2;
		f3 = [l3(ATP)*ActivationStrength*x1^n./(theta_a1^n+x1^n)+l3(ATP)*InhibitionStrength*theta_b1^n./(theta_b1^n+x2^n)-DegradationStrength*x1;...
			l3(ATP)*ActivationStrength*x2^n./(theta_a2^n+x2^n)+l3(ATP)*InhibitionStrength*theta_b2^n./(theta_b2^n+x1^n)-DegradationStrength*x2];
		v3=[x1,x2];
		jac3=jacobian(f3,v3);
			
		for i=0:0.25:5
			for j=0:0.25:5
            matrix_row3=matrix_row3+1;

			X0 = [i,j];
			X = fsolve(fhandle,X0,options);
			x1_ss3=X(1);x2_ss3=X(2);
			x1_ss_3a=round(X(1),3);x2_ss_3a=round(X(2),3);

			sub3=subs(jac3, [x1 x2], [x1_ss3 x2_ss3]);
			eigen3 = eig(sub3);
			eigenvalue_1_3a=sign(eigen3(1));eigenvalue_2_3a=sign(eigen3(2));
            
			sub3b=subs(f3,[x1 x2],[x1_ss3 x2_ss3]);
			
            if (eigenvalue_1_3a < 0) && (eigenvalue_2_3a < 0)
                stability = 1;
            else 
                stability = -1;
            end
			ss_3a(matrix_row3,:) = [DegradationStrength ATP l3(ATP) i j x1_ss_3a x2_ss_3a stability sub3b(1) sub3b(2)];
			end
		end
	col_check=ss_3a(:,9);ss_3a2=ss_3a(col_check < 1e-10,:);
	col_check2=ss_3a2(:,10);ss_3a3=ss_3a2(col_check2 < 1e-10,:);
	col_check3=ss_3a3(:,9);ss_3a4=ss_3a3(col_check3 > -1e-10,:);
	col_check4=ss_3a4(:,10);ss_3a5=ss_3a4(col_check4 > -1e-10,:);
	
	col_k=ss_3a(:,1);ss_3b=ss_3a(col_k==DegradationStrength,:);
	col_stable3=ss_3b(:,8);ss_3c=ss_3b(col_stable3==1,:);ss_3e=unique(ss_3c(:,[1 2 6 7]),'rows');
	col_unstable3=ss_3b(:,8);ss_3d=ss_3b(col_unstable3==-1,:);ss_3f=unique(ss_3d(:,[1 2 6 7]),'rows');
	
	xlim([0 3000]);hold on;ylim([0 4]);grid on;box on;
	plot(ss_3e(:,2),ss_3e(:,4),'b.','MarkerSize',0.7);
	plot(ss_3f(:,2),ss_3f(:,4),'r.','MarkerSize',0.7);
	hx=xlabel('$[ATP] /\mu$M','interpreter','latex');hx.FontSize=fs_labels;hx.FontName=fn;
	hy=ylabel('$x_{2}$ Steady State Position','interpreter','latex');hy.FontSize=fs_labels;hy.FontName=fn;
	ax = gca;ax.FontSize=fs_axis;ax.FontName=fn;ax.TickDir = 'out';
    fig = gcf;fig.Units='centimeters';fig.Position=[0 0 wd ht];fig.PaperUnits='centimeters';fig.PaperSize=[wd ht];
    hold off;
	end
pngFileName5 = sprintf('bif_ATP_k=%.0f.eps',DegradationStrength*100);fullFileName5=fullfile(folder,pngFileName5);print(fig3,fullFileName5,'-depsc');

clock_end = datestr(now,'HH:MM:SS');total_clock = clock_end-clock_start;
fprintf('k=%.2f loop start = %s & end = %s.\n',DegradationStrength,clock_start, clock_end);
end
