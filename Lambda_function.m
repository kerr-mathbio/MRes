fn='Helvetica';fs_labels=10;fs_axis=9;wd=8;ht=7;%wd=10;ht=9;
folder = '/Users/rdk316/Dropbox/PhD/Year 1/Project Rotations/UoB/Matlab & Maple/Varying_ATP/Figures/Function_Plots';
% % folder = '\\phymat.adf.bham.ac.uk\RDK316\PhD\UoB\Figures\Function_Plots';

y=linspace(0,3000,15000);%ATP

d=0.5;e=0.0033;f=-5;g=0.5;%paramter values within lambda function
l= @(y) d*tanh(e*y+f)+g;%lambda([ATP])
%reduce e and scale f accordingly for a shallower function

sigmoid_fig=figure('Name','Sigmoidal Function');clf;
box on;hold on;grid on;
xlim([0 3000]);ylim([0 1]);
plot(y,l(y),'m-');
hx=xlabel('$[ATP]~/\mu$M','interpreter','latex');set(hx,'fontsize',fs_labels);set(hx,'fontname',fn);%x-axis
hy=ylabel('$\lambda$','interpreter','latex');set(hy,'fontsize',fs_labels);set(hy,'fontname',fn);%y-axis
ax = gca;ax.TickDir='out';%ax.XTickLabel=[];ax.YTickLabel=[];
ax.YTick = 0:0.1:1;ax.FontName=fn;ax.FontSize=fs_axis;ax.XTick = 0:500:3000;
set(gcf,'Units','centimeters','Position',[0 0 wd ht],'PaperUnits','centimeters','PaperSize',[wd ht]);
hold off;

pngFileName = 'sigmoidal_function.eps';fullFileName =fullfile(folder, pngFileName);print(sigmoid_fig,fullFileName,'-depsc');
pngFileName2 = 'sigmoidal_function.fig';fullFileName2 = fullfile(folder, pngFileName2);saveas(sigmoid_fig,fullFileName2);

