%%%code in this script is produced with comments explaining what the line of code next to or below it does
%%%code is explained for lambda([ATP]) figure

%pre-setting the font, figure size and fontsizes
fn='Helvetica';fs_labels=10;fs_axis=9;wd=8;ht=7;
%output directory for figure
folder = 'U:\PhD\UoB\Figures';
%ATP concentration range to cover
ATP=linspace(0,3000,15000);
%parameter values for sigmoidal curve
s1=0.5;s2=0.0033;s3=-5;s4=0.5;
%lambda([ATP])
l= @(ATP) s1*tanh(s2*ATP+s3)+s4;

%creating figure
sigmoid_fig=figure('Name','Sigmoidal Function');clf;
%adding grid and box to figure
box on;hold on;grid on;
%x and y axes limits
xlim([0 3000]);ylim([0 1]);
%plotting the sigmoidal curve
plot(ATP,l(ATP),'m-');
%x-axis settings
hx=xlabel('$[ATP]~/\mu$M','interpreter','latex');set(hx,'fontsize',fs_labels);set(hx,'fontname',fn);%x-axis
%y-axis settings
hy=ylabel('$\lambda$','interpreter','latex');set(hy,'fontsize',fs_labels);set(hy,'fontname',fn);%y-axis
%changing x and y axes properties
ax = gca;ax.TickDir='out';ax.YTick = 0:0.1:1;ax.FontName=fn;ax.FontSize=fs_axis;ax.XTick = 0:500:3000;
%setting figure size
set(gcf,'Units','centimeters','Position',[0 0 wd ht],'PaperUnits','centimeters','PaperSize',[wd ht]);
hold off;

%saving produced figure to output directory with specified name and file extenstion
pngFileName = 'sigmoidal_function.eps';fullFileName =fullfile(folder, pngFileName);print(sigmoid_fig,fullFileName,'-depsc');

