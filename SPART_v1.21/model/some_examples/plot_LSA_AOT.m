function plot_LSA_AOT(out_dir)


set(groot,'defaulttextfontsize',15)
set(groot,'defaultAxesfontsize',15)

set(groot,'DefaultTextInterpreter','latex')
set(groot, 'DefaultLegendInterpreter', 'latex')
set(groot,'DefaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultUicontrolFontName', 'Arial')
set(groot, 'defaultUitableFontName', 'Arial')
set(groot, 'defaultAxesFontName', 'Arial')
set(groot, 'defaultTextFontName', 'Arial')
set(groot, 'defaultLegendFontName', 'Arial')
set(groot, 'defaultUipanelFontName', 'Arial')

set(groot,'defaultLineLineWidth',2)
set(groot,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');

cols = [230, 25, 75; 60, 180, 75; 255, 225, 25; 0, 130, 200; 245, 130, 48; 145, 30, 180; 70, 240, 240; 240, 50, 230];
cols = cols/255;

colsx = cols([1:4],:);
%%
directory   = out_dir;
TOC_data    = load([directory,'\TOC_reflectance.dat']);
TOA_data    = load([directory,'\TOA_reflectance.dat']);
wl_sensor   = load([directory,'\wl_spart.dat']);
% params      = readtable([directory,'\pars_and_input.dat']);

% params_array    = table2array(params);


[wl,id_sort]   =   sort(wl_sensor);


if size(TOC_data,1)~=length(wl)
    TOC_data=TOC_data';
end 
if size(TOA_data,1)~=length(wl)
    TOA_data=TOA_data';
end 

TOC_R    =  TOC_data(id_sort,:);
TOA_R    =  TOA_data(id_sort,:);

AOT     =   0:0.5:1.5;
figure(1)
ax1 = axes('Position', [0.1 0.1 0.8 0.7]);
hdlY1 =plot(wl,[TOC_R,TOA_R],'-x');
hold on
axesHandlesToAllLines = findobj(hdlY1, 'Type', 'line');
for j=1:length(AOT)
set(axesHandlesToAllLines(j),'Color',colsx(j,:),'linestyle','--');
set(axesHandlesToAllLines(j+4),'Color', colsx(j,:));
end
key= {'TOC AOT$_{550}$ = 0','TOC AOT$_{550}$ = 0.5','TOC AOT$_{550}$ = 1','TOC AOT$_{550}$ = 1.5',...
    'TOA AOT$_{550}$ = 0','TOA AOT$_{550}$ = 0.5','TOA AOT$_{550}$ = 1','TOA AOT$_{550}$ = 1.5'};
% xlim([400,1100])
ylim([0,0.5])
xlabel('Wavelength (nm)')
ylabel('Reflectance')
oldp = get(gca,'Position');
legend_h = gridLegend(hdlY1,2,key,'location','northoutside','Fontsize',12,'Box','on');
set(gca,'Position',oldp)
set(legend_h,'Position',[0.2,0.82,0.6,0.15])
grid on 
set(gca,'XMinorTick','on','YMinorTick','on')

set(gcf, 'Position',[30,30,1200,700])
plot(wl,0.001,'x','color',[0.6,0.6,0.6])


