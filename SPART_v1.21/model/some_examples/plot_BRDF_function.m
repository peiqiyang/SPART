function plot_BRDF_function(tto,psi,R_TOC,R_TOA,wl,band_i_plot)
% Use: makes BRDF, BFDF and
% bidirectional temperature polar plots from a SPART simulations
% of directional data.
if size(tto,1)~=1
    obs_zenith                  = tto';
end
if size(psi,1)~=1
    obs_azimuth                  = psi';
end
wl_i                        =  wl(band_i_plot);

obs_azimuth                 = obs_azimuth-360*(obs_azimuth>180);
% obs_zenith_i                = 0:90;
obs_zenith_i                = 0:70;
obs_azimuth_i               = -180:1:180;
[Obs_Zenith_i,Obs_Azimuth_i]= meshgrid(obs_zenith_i,obs_azimuth_i);

x                           = obs_zenith  *pi/180 .* cos(obs_azimuth  *pi/180+pi/2);
y                           = obs_zenith  *pi/180 .* sin(obs_azimuth  *pi/180+pi/2);

X_i                         = Obs_Zenith_i*pi/180 .* cos(Obs_Azimuth_i*pi/180+pi/2);
Y_i                         = Obs_Zenith_i*pi/180 .* sin(Obs_Azimuth_i*pi/180+pi/2);

R_TOC_i                     = griddata(x,y,R_TOC(band_i_plot,:),X_i,Y_i,'v4');
R_TOA_i                     = griddata(x,y,R_TOA(band_i_plot,:),X_i,Y_i,'v4');


%%
xli = .5*pi*[0 -1  -.05 0.8];
yli = .5*pi*[0.9  0 -0.9 0];

ax1 = axes('Position', [0.1,0.6,0.24,0.25]);
z = pcolor(X_i,Y_i,R_TOC_i); hold on
set(z,'LineStyle','none')
for k = 1:4
    plot(20*k/180*pi.*cos((-1:.01:1)*pi),20*k/180*pi.*sin((-1:.01:1)*pi),'Color',[1 1 1])
    text(20*k/180*pi-.2,.2,num2str(20*k),'FontSize',14,'Color',[0 0 0]);
    text(xli(k),yli(k),num2str(90*(k-1)),'FontSize',14,'Color','k','FontAngle','italic');
end
    text(-0.8,1.8,['TOC R (',num2str(wl_i),'nm)'],'FontSize',14)
        colormap jet
    h = colorbar;
axis off
set(h,'position',get(h,'position')+[0.15,0,0,0])


ax2 = axes('Position', [0.54,0.6,0.24,0.25]);
z = pcolor(X_i,Y_i,R_TOA_i); hold on
set(z,'LineStyle','none')
for k = 1:4
    plot(20*k/180*pi.*cos((-1:.01:1)*pi),20*k/180*pi.*sin((-1:.01:1)*pi),'Color',[1 1 1])
    text(20*k/180*pi-.2,.2,num2str(20*k),'FontSize',14,'Color',[0 0 0]);
    text(xli(k),yli(k),num2str(90*(k-1)),'FontSize',14,'Color','k','FontAngle','italic');
end
 text(-0.8,1.8,['TOA R (',num2str(wl_i),'nm)'],'FontSize',14)
 colormap jet
 h = colorbar;
axis off
set(h,'position',get(h,'position')+[0.15,0,0,0])
%% solar plane 
psi1     = obs_azimuth(obs_azimuth ==0);
psi2     = obs_azimuth(obs_azimuth ==180);

tto1     =   obs_zenith(obs_azimuth ==0);
tto2     = - obs_zenith(obs_azimuth ==180);



R_TOCi1   = R_TOC(band_i_plot,obs_azimuth ==0);
R_TOCi2   = R_TOC(band_i_plot,obs_azimuth ==180);

R_TOAi1   = R_TOA(band_i_plot,obs_azimuth ==0);
R_TOAi2   = R_TOA(band_i_plot,obs_azimuth ==180);


% psi  = [psi1,psi2];
tto  = [tto1,tto2];
R_TOCi  = [R_TOCi1,R_TOCi2];
R_TOAi  = [R_TOAi1,R_TOAi2];

[tto,ord]   =   sort(tto);
R_TOCi      =   R_TOCi(ord);
R_TOAi      =  R_TOAi(ord);


ax3 = axes('Position', [0.1,0.1,0.8,0.3]);
plot(tto,R_TOCi,'x--','linewidth',2,'color',[0.7,0.7,0.7])
hold on
plot(tto,R_TOAi,'x--','linewidth',2,'color',[0.7,0.7,0.7])
legend('TOC','TOA','location','north')
xlim([-70,70])
title('Reflectance at principal plane')
xlabel('VZA')
ylabel('Reflectance')
grid on 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gcf, 'Units','centimeters','Position',[3,3,20,18])

