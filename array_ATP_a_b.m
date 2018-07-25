%%%code in this script is produced with comments explaining what the line of code next to or below it does.
%%%code is explained for all phase portrait array figures

%pre-setting the font, figure size and fontsizes
fn='Helvetica';wd=16;ht=14;fs_labels=10;fs_axis=8;
%parameter values for ODEs that are being fixed
DegradationStrength=1;n=4;theta_a1=0.5;theta_a2=0.5;theta_b1=0.5;theta_b2=0.5;
%parameter values for sigmoidal curve
s1=0.5;s2=0.0033;s3=-5;s4=0.5;
%fsolve tolerances
options = optimoptions('fsolve','FunctionTolerance',1e-11,'OptimalityTolerance',1e-11,'Display','off');
%output directory for figures
folder = 'U:\PhD\UoB\Figures\Num_SS_Arrays';
%setting the initial value of the row in matrix ss_matrix as zero
matrix_row=0;
%pre-setting matrix size to speed up computations
ss_matrix=zeros(95000,10);

%pairs of a and b values to be used for arrays
%%% (a,b)=(1,1),(2,1),(1,2),(1,0.25),(0.25,1)
for InhibitionStrength=1 %set b value(s)
    for ActivationStrength=2 %set a value(s)
		%starting time of loop
        clock_start = datestr(now,'HH:MM:SS');
		%display where the computation is up to in command window - good for long computations to see where you are up to
		fprintf('Start of ActivationStrength=%.2f, InhibitionStrength=%.2f computation at %s.\n',ActivationStrength,InhibitionStrength,datestr(now,'HH:MM:SS'));
		%creating figure
		array_fig=figure;
		%setting figure size
        set(gcf,'Units','centimeters','Position',[0 0 wd ht],'PaperUnits','centimeters','PaperSize',[wd ht]);
		%setting value of sub-figure position to zero
		sub_figure_position=0;
		for ATP = [1000 1200 1300 1350 1400 1450 1500 1550 1600 1650 1700 1750 1800 1900 2000 2700]%ATP values to scan through
			%increasing the value by 1 in each new loop - sets sub-figure position
			sub_figure_position=sub_figure_position+1;
			%lambda([ATP])
			l= @(ATP) s1*tanh(s2*ATP+s3)+s4;
			%display where the computation is up to in command window - good for long computations to see where you are up to
			fprintf('Running ActivationStrength=%.2f, InhibitionStrength=%.2f and [ATP]=%d computation at %s.\n',ActivationStrength,InhibitionStrength,ATP,datestr(now,'HH:MM:SS'));
			%setting number of sub-figures and sub-figure position
			subplot(4,4,sub_figure_position);
            
			%symbolic variables for protein levels
            syms x1 x2;
			%ODEs
            f1 = [l(ATP)*ActivationStrength*x1^n./(theta_a1^n+x1^n)+l(ATP)*InhibitionStrength*theta_b1^n./(theta_b1^n+x2^n)-DegradationStrength*x1;...
                l(ATP)*ActivationStrength*x2^n./(theta_a2^n+x2^n)+l(ATP)*InhibitionStrength*theta_b2^n./(theta_b2^n+x1^n)-DegradationStrength*x2];
            % variables for jacobian matrix
			v=[x1,x2];
			%calculating jacbian with respect to variables x1 & x2
            jac=jacobian(f1,v);
            %function to use by fsolve 
            fhandle=@(X)vary_ActivationStrength_InhibitionStrength_ATP(X,ActivationStrength,InhibitionStrength,ATP);
            
			for i=0:0.25:5%initial conditions on x axis
				for j=0:0.25:5%initial conditions on y axis
					%increasing value of the row in matrix ss_matrix by 1 in each loop
					matrix_row=matrix_row+1;
					
                    %initial conditions pairing
					X0 = [i,j];
					%using fsolve
					X = fsolve(fhandle,X0,options);
					%calculated steady state values
					x1_ss=X(1);x2_ss=X(2);
					%rounded calculated steady state values
					x1_ss_b=round(X(1),3);x2_ss_b=round(X(2),3);
					
					%subs. in steady state values to jacobian
					sub=subs(jac, [x1 x2], [x1_ss x2_ss]); 
					%calc eigenvlaues
					eigen = eig(sub); 
					%calculate the sign of each eigenvalue
					eigenvalue_1_2=sign(eigen(1));eigenvalue_2_2=sign(eigen(2));
					%substituting the calculated steady state back into ODEs
					sub2=subs(f1,[x1 x2],[x1_ss x2_ss]);
					
					%testing if the steady state is stable or unstable
					if (eigenvalue_1_2 < 0) && (eigenvalue_2_2 < 0)
						stability = 1;
					else 
						stability = -1;
					end
					
					%matrix of a, b, ATP, ics, ss positions, the stability and subs2
					ss_matrix(matrix_row,:) = [ActivationStrength InhibitionStrength ATP i j x1_ss_b x2_ss_b stability sub2(1) sub2(2)];
					
					%checking if the steady state when sub into ODEs is within a tolerance of plus-minus 10^-10
					col_check=ss_matrix(:,9);ss_b=ss_matrix(col_check < 1e-10,:);
					col_check2=ss_b(:,10);ss_c=ss_b(col_check2 < 1e-10,:);
					col_check3=ss_c(:,9);ss_d=ss_c(col_check3 > -1e-10,:);
					col_check4=ss_d(:,10);ss_1=ss_d(col_check4 > -1e-10,:);
					
					%selecting the steady states for the current values of a, b and ATP
					col_a1=ss_1(:,1);%selecting a column in ss_1 matrix
					ss_1b = ss_1(col_a1 == ActivationStrength,:);%new sub-matrix ss_1b is a submatrix of ss_1 with value a=a in the first column
					col_b1=ss_1b(:,2);%selecting a column in ss_1b matrix
					ss_1c = ss_1b(col_b1 == InhibitionStrength,:);%new sub-matrix ss_1c is a submatrix of ss_1b with value b=b in the second column
					col_y1 = ss_1c(:,3);%selecting a column in ss_1c matrix
					ss_1d = ss_1c(col_y1 == ATP,:);%new matrix ss_1d is a submatrix of ss_1c with value y=y in the third column
					
					%extracting stable and unstable sub-matrices
					stab_col1=ss_1d(:,8);stable_ss1 = ss_1d(stab_col1==1,:);unique_stable_ss1 = unique(stable_ss1(:,[1 2 3 6 7]),'rows');
					unstab_col1=ss_1d(:,8);unstable_ss1 = ss_1d(stab_col1==-1,:);unique_unstable_ss1 = unique(unstable_ss1(:,[1 2 3 6 7]),'rows');
					
					%plotting sub-figre 
					grid on;hold on;box on; %axes limits & box around figure
					plot(unique_unstable_ss1(:,4),unique_unstable_ss1(:,5),'ro','MarkerSize',1.8,'MarkerFaceColor','r');%plotting stable steady states
					plot(unique_stable_ss1(:,4),unique_stable_ss1(:,5),'bo','MarkerSize',1.8,'MarkerFaceColor','b');%plotting unstable steady states
                    ax = gca;ax.XTick = 0:1:4;ax.YTick = 0:1:4;ax.YLim = [0 4];ax.XLim = [0 4];%changing x and y axes properties
                    ax.FontSize=fs_axis;ax.FontName=fn;ax.TickDir = 'out';%changing fontsize, fontname and tick direction on current sub-figure
                    hx=xlabel('{  }');hx.Interpreter='latex';hx.FontSize=fs_labels;hx.FontName=fn;%no label on x-axis
                    hy=ylabel('{  }');hy.Interpreter='latex';hy.FontSize=fs_labels;hy.FontName=fn;%no label on y-axis
                    hold off;
				end
			end
		end
		%saving produced figure to output directory with specified name and file extenstion
	epsFileName = sprintf('Array_ATP_a=%.0f_b=%.0f.eps',ActivationStrength*100,InhibitionStrength*100);fullFileName =fullfile(folder, epsFileName);print(array_fig,fullFileName,'-depsc');
    
	%end time of loop
    clock_end = datestr(now,'HH:MM:SS');
	%displays start and end time of loop computation in command window
	fprintf('ActivationStrength=%.2f, InhibitionStrength=%.2f loop start = %s & end = %s.\n',ActivationStrength,InhibitionStrength,clock_start, clock_end);
    end
end