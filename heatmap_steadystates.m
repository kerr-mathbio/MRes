%%%code in this script is produced with comments explaining what the line of code next to or below it does
%%%code is explained for all heatmap figures

%pre-setting the font, figure size and fontsizes
fn='Helvetica';wd=10;ht=9;fs_labels=10;fs_axis=9;
%parameter values for ODEs that are being fixed
DegradationStrength=1;n=4;theta_a1=0.5;theta_a2=0.5;theta_b1=0.5;theta_b2=0.5;
%parameter values for sigmoidal curve
s1=0.5;s2=0.0033;s3=-5;s4=0.5;
%fsolve tolerances
options = optimoptions('fsolve','FunctionTolerance',1e-11,'OptimalityTolerance',1e-11,'Display','off'); 
%output directory for figures
folder = 'U:\PhD\UoB\Figures\Heatmaps';
%setting the initial row in matrix 'ss_matrix' as zero
matrix_row=0;%q=0;
%pre-setting matrix size to speed up computations
ss_matrix=zeros(1200000,11);
%order of ATP concetrations for y-axis -- see use further down
NO = {'2700','2000','1900','1800','1750','1700','1650','1600','1550','1500','1450','1400','1350','1300','1200','1000'};

for ActivationStrength=[0,0.5,1,1.5,2,3]%set a values
	%starting time of loop
	clock_start = datestr(now,'HH:MM:SS');
	%display where the computation is up to in command window - good for long computations to see where you are up to
	fprintf('ActivationStrength=%.2f heatmap start time = %s\n.',ActivationStrength,datestr(now,'HH:MM:SS'));
	for InhibitionStrength=0:0.25:3%set b values
        for y = [1000 1200 1300 1350 1400 1450 1500 1550 1600 1650 1700 1750 1800 1900 2000 2700]%ATP values to scan through
			%lambda
			l= @(ATP) s1*tanh(s2*ATP+s3)+s4;
			%display where the computation is up to in command window - good for long computations to see where you are up to
			fprintf('Running heatmap with ActivationStrength=%.2f, InhibitionStrength=%.2f and [ATP]=%d at %s.\n',ActivationStrength,InhibitionStrength,y,datestr(now,'HH:MM:SS'));
			%fsolve function
			fhandle=@(X)bif_a_b_y(X,ActivationStrength,InhibitionStrength,y);
			%symbolic variables for protein levels
			syms x1 x2;
			%ODEs
			f1 = [l(ATP)*ActivationStrength*x1^n./(theta_a1^n+x1^n)+l(ATP)*InhibitionStrength*theta_b1^n./(theta_b1^n+x2^n)-DegradationStrength*x1;...
				l(ATP)*ActivationStrength*x2^n./(theta_a2^n+x2^n)+l(ATP)*InhibitionStrength*theta_b2^n./(theta_b2^n+x1^n)-DegradationStrength*x2];
			%variables for jacobian matrix
			v=[x1,x2];
			%calculating jacbian with respect to variables x1 & x2
			jac=jacobian(f1,v);
					
            for i=0:0.25:5%initial conditions on x axis
				for j=0:0.25:5%initial conditions on y axis
					%moving row in matrix down by 1 in each loop
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
					sub=subs(jac, [x1 x2], [x1_ss x2_ss]); %subs. in ss values from original ics
					%calc eigenvlaues
					eigen = eig(sub); %calc eigenvlaues of matrix 'sub'
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
					
					%rounding lambda
					lambda = round(l(y),3);
					%matrix of a, b, ATP, ics, ss positions, the stability, lambda and subs2
					ss_matrix(matrix_row,:) = [ActivationStrength InhibitionStrength y i j x1_ss_b x2_ss_b stability lambda sub2(1) sub2(2)];%matrix of value of a, b, ics, ss positions for x1, x2 and the stability of the ss point
				end
            end
        end
	end
%checking if the steady state when sub into ODEs is within a tolerance of plus-minus 10^-10
col_check=ss_matrix(:,10);ss_b=ss_matrix(col_check < 1e-10,:);
col_check2=ss_b(:,11);ss_c=ss_b(col_check2 < 1e-10,:);
col_check3=ss_c(:,10);ss_d=ss_c(col_check3 > -1e-10,:);
col_check4=ss_d(:,11);ss_1=ss_d(col_check4 > -1e-10,:);
					
%selecting a column in ss_1 matrix
stab_col3=ss_1(:,8);
%new sub-matrix ss_8 is a submatrix of ss_1 with stable steady states
ss_8=ss_1(stab_col3 == 1,:);
%extracting unique stable steady states
ss_9 = unique(ss_8(:,[1 2 3 6 7 9]),'rows');
%selecting a column in ss_9 matrix
col_a4 = ss_9(:,1);
%new sub-matrix stable_ss3 is a submatrix of ss_9 with a=a
stable_ss3 = ss_9(col_a4 == ActivationStrength,:);

%table of a, b, ATP, stable steady state positions and lambda value
T2 = array2table(stable_ss3,'VariableNames',{'a','b','ATP','stable_ss_position_x1','stable_ss_position_x2','lambda'});

%creating figure
fig_heatmap = figure('Name','Heatmap');
%plotting heatmap from table T2' with b on the x-axis and ATP values on the y-axis
h = heatmap(T2,'b','ATP','Title','');
%colour scheme to use -- could use spring or summer or parula
h.Colormap = cool;
%data colour and label if missing
h.MissingDataColor = [0.8 0.8 0.8];h.MissingDataLabel = 'No data';
%method to use for displaying data
h.ColorMethod = 'count';
%color bar visible or not in figures 
h.ColorbarVisible = 'off';
%fontname, colour and fontsize
h.FontColor = 'black';h.FontName = fn;h.FontSize = fs_labels;
%colour limits to keep consistency between figures 
h.ColorLimits=[1 4];

%re-arranging y-axis ATP values
h.SourceTable.ATP = categorical(h.SourceTable.ATP);
neworder = NO;
h.SourceTable.ATP = reordercats(h.SourceTable.ATP,neworder);

%removing y-axis label
h.YLabel=' ';
%removing x-axis label
h.XLabel=' ';
%axis fontname and fontsize
ax = gca;ax.FontSize=fs_axis;ax.FontName=fn;
%figure size
fig = gcf;fig.Units='centimeters';fig.Position=[0 0 wd ht];fig.PaperUnits='centimeters';fig.PaperSize=[wd ht];

%saving produced figure to output directory with specified name and file extenstion
pngFileName3 = sprintf('Heatmap_ATP_a=%.0f.eps',ActivationStrength*100);fullFileName3 = fullfile(folder, pngFileName3);print(fig_heatmap,fullFileName3,'-depsc');

%end time of loop
clock_end = datestr(now,'HH:MM:SS');
%displays start and end time of loop computation in command window
fprintf('a=%.2f heatmap start = %s & end = %s.\n',ActivationStrength,clock_start,clock_end);
end
   