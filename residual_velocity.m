clear all
close all
clc

% addpath /Users/janneke/Documents/openearthtools2/openearthtools/trunk/matlab 
% oetsettings

run = 'Hydro_RZcomp';
dx=10;
dt=15;

workdir = ['/Volumes/SEAGATE_EXP/2018-2019/Marid/',num2str(run),'/'];
figdir = '/Volumes/SEAGATE_EXP/2018-2019/Marid/Figures/';

RZ =1;
RR =0;
%%
GridFile= fullfile(['/Volumes/SEAGATE_EXP/2018-2019/Marid/',num2str(run),'/grid_10m_v3.grd']);
DepthFile = fullfile(['/Volumes/SEAGATE_EXP/2018-2019/Marid/',num2str(run),'/dep_2000_2.dep']);

GRID = wlgrid('read',GridFile);
depth = wldep('read',DepthFile,GRID);

xa = 46;
xb = 533;
%%
trimfile1=fullfile(workdir,'trim-bas_sed.dat');
trim1 = vs_use(trimfile1);
m1 = qpfopen(trim1.DatExt);
l1 = qpread(m1);

%%
xZ = 2;
xN = length(depth(:,2))-1;

%%
time=vs_time(trim1);
Morfotime=time.t/60;
ttime=Morfotime*60;
Morfotime3=time.datenum;

%%
U1 = vs_get(trim1,'map-series',{1:length(Morfotime)},'U1',{2 0 0});
W1 = vs_get(trim1,'map-series',{1:length(Morfotime)},'WPHY',{2 0 0});
x = vs_get(trim1,'map-const',{1},'XCOR',{2 0});
y = vs_get(trim1,'map-const',{1},'YCOR',{2 0});
dz = vs_get(trim1,'map-const',{1},'THICK',{0});
wl1 = vs_get(trim1,'map-series',{1:length(Morfotime)},'S1',{2 0});
Qs1 = vs_get(trim1,'map-sed-series',{1:length(Morfotime)},'SSUU',{2 0 1});
Qb1 = vs_get(trim1,'map-sed-series',{1:length(Morfotime)},'SBUU',{2 0 1});
Zb1 = vs_get(trim1,'map-sed-series',{1:length(Morfotime)},'DPS',{2 0});
%%
u=cell2mat(U1);
w=cell2mat(W1);
wl=cell2mat(wl1);
Qs=cell2mat(Qs1);
Qb=cell2mat(Qb1);
Zb=cell2mat(Zb1);
%%
depth = Zb+wl;
dz2=cumsum(dz);
for aa=1:length(depth(:,1)); 
    for bb=1:length(depth(1,:))-1;
        for cc=1:length(dz);
            z(aa,bb,cc)=dz2(cc)*-depth(aa,bb);
        end
    end
end

distance = sqrt((x(1:end-1)-x(1)).^2+(y(1:end-1)-y(1)).^2);

%%
if RZ==1;
    ta = 6;
    tb = 451;
    no = 9;
elseif RR==1;
    ta = 9;
    tb = 454;
    no = 9;
end

%%
figure;
plot(u(:,245,1))
hold on
scatter(ta,u(ta,245,1))
scatter(tb,u(tb,245,1))

%%
for a = 1:length(u(1,:,1))-1;
    for b = 1:length(u(1,1,:));
%         Ures(a,b) = (1/(length(ta:tb)*dt*60))*sum(u(ta:tb,a,b)*dt*60);
        Ures(a,b) = nanmean(u(ta:tb,a,b));
%         Wres(a,b) = (1/(length(ta:tb)*dt*60))*sum(w(ta:tb,a,b)*dt*60);
        Wres(a,b) = nanmean(w(ta:tb,a,b));
    end
end


for jj=1:length(Qb(1,:))-1;
    Qbres(jj)=nanmean(Qb(ta:tb,jj));
    Qsres(jj)=nanmean(Qs(ta:tb,jj));
end

Qtres=Qbres+Qsres;

dQb = gradient(Qbres,distance);
dQs = gradient(Qsres,distance);
dQt = gradient(Qtres,distance);
%%
for c = 1:length(u(1,1,:));
    Ures2(:,c) = Ures(:,c)-nanmean(Ures(:,c));
end

%%
for ii=1:50;
    xx(:,ii)=distance;
end
zz=squeeze(z(15,:,:));

%%
save([run,'_resflow.mat'],'u','Morfotime','Ures','Wres','Ures2','Qbres','Qsres','Qtres','dQb','dQs','dQt','Zb','distance','xx','zz');
%%
levelsU=[0:5e-3:5e-2];

figure;
contourf(xx/1E3,zz,Ures,300,'LineStyle','none')
hold on
[C,h]=contour(xx/1E3,zz,Ures,levelsU,'ShowText','on');
h.LineColor='k';
clabel(C,h,'FontSize',15,'FontName','Times New Roman')
caxis([0 5e-2])
set(gca,'FontSize',22,'FontName','Times New Roman')
b = colorbar;
b.Label.String = '\langleu\rangle [m s^{-1}]';
xlabel('x [km]')
ylabel('z [m]')
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1400,400]);
set(gcf,'PaperPositionMode','auto')
eval(['print -depsc -tiff ' ,[figdir,run,'_Ures.eps']])

figure;
contourf(xx(xa:xb,:)/1E3,zz(xa:xb,:),Ures(xa:xb,:),300,'LineStyle','none')
hold on
[C,h]=contour(xx(xa:xb,:)/1E3,zz(xa:xb,:),Ures(xa:xb,:),levelsU,'ShowText','on');
h.LineColor='k';
clabel(C,h,'FontSize',15,'FontName','Times New Roman')
caxis([0 5e-2])
set(gca,'FontSize',22,'FontName','Times New Roman')
b = colorbar;
b.Label.String = '\langleu\rangle [m s^{-1}]';
xlabel('x [km]')
ylabel('z [m]')
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1400,400]);
set(gcf,'PaperPositionMode','auto')
eval(['print -depsc -tiff ' ,[figdir,run,'_Ures_zoom.eps']])

%%
levelsU2=[-1e-2:2e-3:1e-2];

figure;
contourf(xx(xa:xb,:)/1E3,zz(xa:xb,:),Ures2(xa:xb,:),300,'LineStyle','none')
hold on
[C,h]=contour(xx(xa:xb,:)/1E3,zz(xa:xb,:),Ures2(xa:xb,:),levelsU2,'ShowText','on');
h.LineColor='k';
clabel(C,h,'FontSize',15,'FontName','Times New Roman')
caxis([-1e-2 1e-2])
set(gca,'FontSize',22,'FontName','Times New Roman')
b = colorbar;
b.Label.String = '\langleu\rangle [m s^{-1}]';
xlabel('x [km]')
ylabel('z [m]')
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1400,400]);
set(gcf,'PaperPositionMode','auto')
eval(['print -depsc -tiff ' ,[figdir,run,'_Ures2_zoom.eps']])
%%
levelsW=[-1e-3:5e-4:1e-3];

figure;
contourf(xx/1E3,zz,Wres,300,'LineStyle','none')
hold on
[C,h]=contour(xx/1E3,zz,Wres,levelsW,'ShowText','on');
h.LineColor='k';
clabel(C,h,'FontSize',15,'FontName','Times New Roman')
caxis([-1e-3 1e-3])
set(gca,'FontSize',22,'FontName','Times New Roman')
b = colorbar;
b.Label.String = '\langlew\rangle [m s^{-1}]';
xlabel('x [km]')
ylabel('z [m]')
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1400,400]);
set(gcf,'PaperPositionMode','auto')
eval(['print -depsc -tiff ' ,[figdir,run,'_Wres.eps']])

figure;
contourf(xx(xa:xb,:)/1E3,zz(xa:xb,:),Wres(xa:xb,:),300,'LineStyle','none')
hold on
[C,h]=contour(xx(xa:xb,:)/1E3,zz(xa:xb,:),Wres(xa:xb,:),levelsW,'ShowText','on');
h.LineColor='k';
clabel(C,h,'FontSize',15,'FontName','Times New Roman')
caxis([-1e-3 1e-3])
set(gca,'FontSize',22,'FontName','Times New Roman')
b = colorbar;
b.Label.String = '\langlew\rangle [m s^{-1}]';
xlabel('x [km]')
ylabel('z [m]')
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1400,400]);
set(gcf,'PaperPositionMode','auto')
eval(['print -depsc -tiff ' ,[figdir,run,'_Wres_zoom.eps']])

%%
fig=figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(distance(xa:xb)/1E3,Qbres(xa:xb),'-r','LineWidth',2)
hold on
plot(distance(xa:xb)/1E3,Qsres(xa:xb),'Color',[0 0.7 0],'LineStyle','-','LineWidth',2)
plot(distance(xa:xb)/1E3,Qtres(xa:xb),'-b','LineWidth',2)
set(gca,'FontSize',22,'FontName','Times New Roman')
h=legend('Bed load','Suspended load','Total load','location','southwest')
set(h.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));
ylabel('\langleQ\rangle [m^{3}m^{-1}s^{-1}]')
xlim([21.5 25.5])
ylim([0 1.5E-6])
yyaxis right
plot(distance(xa:xb)/1E3,-Zb(1,xa:xb),':k','LineWidth',2,'HandleVisibility','off')

set(gca,'FontSize',22,'FontName','Times New Roman')
ylabel('z_b [m]')
xlabel('x [km]')

set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1400,400]);
set(gcf,'PaperPositionMode','auto')
eval(['print -depsc -tiff ' ,[figdir,run,'_Qres_zoom.eps']])

%%
fig=figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(distance(xa:xb)/1E3,smooth(dQb(xa:xb)),'-r','LineWidth',2)
hold on
plot(distance(xa:xb)/1E3,smooth(dQs(xa:xb)),'Color',[0 0.7 0],'LineStyle','-','LineWidth',2)
plot(distance(xa:xb)/1E3,smooth(dQt(xa:xb)),'-b','LineWidth',2)
set(gca,'FontSize',22,'FontName','Times New Roman')
h=legend('Bed load','Suspended load','Total load','location','southwest')
set(h.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));
ylabel('d\langleQ\rangle/dx [m s^{-1}]')
xlim([21.5 25.5])
ylim([-15E-9 5E-9])
yyaxis right
plot(distance(xa:xb)/1E3,-Zb(1,xa:xb),':k','LineWidth',2,'HandleVisibility','off')

set(gca,'FontSize',22,'FontName','Times New Roman')
ylabel('z_b [m]')
xlabel('x [km]')

set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1400,400]);
set(gcf,'PaperPositionMode','auto')
eval(['print -depsc -tiff ' ,[figdir,run,'_dQ_zoom.eps']])