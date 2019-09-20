clear all
close all
clc

addpath /Users/janneke/Documents/openearthtools2/openearthtools/trunk/matlab 
oetsettings

xa = 46;
xb = 533;
%%
GridFile= fullfile(['/Volumes/SEAGATE_EXP/2018-2019/Marid/Timeseries/RZtime_hc_MF50/grid_10m_v3.grd']);
GRID = wlgrid('read',GridFile);

DepthFile = fullfile(['/Volumes/SEAGATE_EXP/2018-2019/Marid/Timeseries/RRtime_hc_MF50/dep_2000_2.dep']);
depth = wldep('read',DepthFile,GRID);
% 
% X = GRID.X(:,2);
% Y = GRID.Y(:,2);
% Z = -depth(1:end-1,2);
% save('Afgetopt.mat','X','Y','Z')

%% data 2010
DepthFile10 = fullfile('/Users/janneke/Documents/2018-2019/Marid/','dep_2010_2.dep');
depth10 = wldep('read',DepthFile10,GRID);

%%
workdir1 = '/Volumes/SEAGATE_EXP/2018-2019/Marid/Timeseries/RRtime_hc_MF50/';
trimfile1=fullfile(workdir1,'trim-bas_sed.dat');
trim1 = vs_use(trimfile1);
m1 = qpfopen(trim1.DatExt);
l1 = qpread(m1);

[succes,Zb1] = qp_getdata(trim1,'bed level in water level points','griddata',0,0,0);
% [succes,wl1] = qp_getdata(trim1,'water level','griddata',0,0,0);
%%
workdir2 = '/Volumes/SEAGATE_EXP/2018-2019/Marid/Timeseries/RZtime_hc_MF50/';
trimfile2=fullfile(workdir2,'trim-bas_sed.dat');
trim2 = vs_use(trimfile2);
m2 = qpfopen(trim2.DatExt);
l2 = qpread(m2);

[succes,Zb2] = qp_getdata(trim2,'bed level in water level points','griddata',0,0,0);
% [succes,wl2] = qp_getdata(trim2,'water level','griddata',0,0,0);

%%
workdir3 = '/Volumes/SEAGATE_EXP/2018-2019/Marid/Timeseries/UZtime_hc_MF50/';
trimfile3=fullfile(workdir3,'trim-bas_sed.dat');
trim3 = vs_use(trimfile3);
m3 = qpfopen(trim3.DatExt);
l3 = qpread(m3);

[succes,Zb3] = qp_getdata(trim3,'bed level in water level points','griddata',0,0,0);
% [succes,wl3] = qp_getdata(trim3,'water level','griddata',0,0,0);

%%
workdir4 = '/Volumes/SEAGATE_EXP/2018-2019/Marid/Timeseries/ZZtime_hc_MF50/';
trimfile4=fullfile(workdir4,'trim-bas_sed.dat');
trim4 = vs_use(trimfile4);
m4 = qpfopen(trim4.DatExt);
l4 = qpread(m4);

[succes,Zb4] = qp_getdata(trim4,'bed level in water level points','griddata',0,0,0);
% [succes,wl4] = qp_getdata(trim3,'water level','griddata',0,0,0);
%%
time=vs_time(trim1);
temp=vs_get(trim1,'map-infsed-serie','MORFT');
Morfotime=cell2mat(temp);

dt = Morfotime(end)*3600*24-Morfotime(1)*3600*24;

%%
time2=vs_time(trim2);
temp2=vs_get(trim2,'map-infsed-serie','MORFT');
Morfotime2=cell2mat(temp2);

dt2 = Morfotime2(end)*3600*24-Morfotime2(1)*3600*24;
%%
x1 = Zb1.X(:,2);
x2 = Zb2.X(:,2);
x3 = Zb3.X(:,2);
x4 = Zb4.X(:,2);

y1 = Zb1.Y(:,2);
y2 = Zb2.Y(:,2);
y3 = Zb3.Y(:,2);
y4 = Zb4.Y(:,2);

dist1 = sqrt((x1-x1(2)).^2+(y1-y1(2)).^2)/1E3;
dist2 = sqrt((x2-x2(2)).^2+(y2-y2(2)).^2)/1E3;
dist3 = sqrt((x3-x3(2)).^2+(y3-y3(2)).^2)/1E3;
dist4 = sqrt((x4-x4(2)).^2+(y4-y4(2)).^2)/1E3;

bed1 = squeeze(Zb1.Val(:,:,2));
bed2 = squeeze(Zb2.Val(:,:,2));
bed3 = squeeze(Zb3.Val(:,:,2));
bed4 = squeeze(Zb4.Val(:,:,2));

% zeta1 = squeeze(wl1.Val(:,:,2));
% zeta2 = squeeze(wl2.Val(:,:,2));

%%
bed1(end,273:285)=bed1(end,273)+(bed1(end,285)-bed1(end,273))/(length(273:285)*10)*(0:10:10*length(273:285-1));
bed3(end,273:292)=bed3(end,273)+(bed3(end,292)-bed3(end,273))/(length(273:292)*10)*(0:10:10*length(273:292-1));
bed4(end,273:291)=bed4(end,273)+(bed4(end,291)-bed4(end,273))/(length(273:291)*10)*(0:10:10*length(273:291-1));

figure;
plot(bed1(end,:))
hold on
plot(bed3(end,:))
plot(bed4(end,:))
%%
figure;
plot(dist1(xa:xb),bed1(1,xa:xb),':k','LineWidth',2)
hold on
% plot(dist1(xa:xb),bed1(floor(end/4),xa:xb),'r','LineWidth',2)
% plot(dist1(xa:xb),bed1(floor(end/2),xa:xb),'g','LineWidth',2)
% plot(dist1(xa:xb),bed1(3*floor(end/4),xa:xb),'m','LineWidth',2)
plot(dist1(xa:xb),bed1(end,xa:xb),'b','LineWidth',2)
plot(dist2(xa:xb),bed2(end,xa:xb),'r','LineWidth',2)
plot(dist3(xa:xb),bed3(end,xa:xb),'Color',[0 0.7 0],'LineWidth',2)
plot(dist4(xa:xb),bed4(end,xa:xb),'m','LineWidth',2)
plot(dist2(xa:xb),-depth10(xa:xb,2),'k','LineWidth',2)
% plot(dist2(xa:xb),bed1(end-2,xa:xb),'b','LineWidth',2)
xlim([20.5 25.5])
set(gca,'FontSize',22,'FontName','Times New Roman')
h=legend('bed level 1999','RR','RZ','UZ','ZZ','data 2010','location','southwest')
set(h.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));
xlim([21 24])
xlabel('x [km]')
ylabel('z_b [m]')
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1800,400]);
set(gcf,'PaperPositionMode','auto')

%%
fig = figure;
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(dist2(xa:xb),bed2(end-2,xa:xb)+depth10(xa:xb,2)','-b','LineWidth',2)
hold on
plot(dist1(xa:xb),bed1(end-1,xa:xb)+depth10(xa:xb,2)','-r','LineWidth',1)
set(gca,'FontSize',22,'FontName','Times New Roman')
ylabel('\Delta z_b [m]')
xlabel('x [km]')
xlim([20.5 25.5])
ylim([-2.5 2.5])


yyaxis right
plot(dist2(xa:xb),-depth10(xa:xb,2)',':k','LineWidth',2)
legend('(MF25) - (data)','(MF50) - (data)','data 2010')
ylim([-25 -20])
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1800,400]);
set(gcf,'PaperPositionMode','auto')
%%
RMSD_RR = RootMeanSquareError(bed1(end,xa:xb),-depth10(xa:xb,2)',length(xa:xb));
RMSD_RZ = RootMeanSquareError(bed2(end,xa:xb),-depth10(xa:xb,2)',length(xa:xb));
RMSD_UZ = RootMeanSquareError(bed3(end,xa:xb),-depth10(xa:xb,2)',length(xa:xb));
RMSD_ZZ = RootMeanSquareError(bed4(end,xa:xb),-depth10(xa:xb,2)',length(xa:xb));

%%
[pks1, locs1] = findpeaks(bed1(end,xa:xb),dist1(xa:xb),'MinPeakDistance',0.25);
[dal1, locsd1] = findpeaks(-bed1(end,xa:xb),dist1(xa:xb),'MinPeakDistance',0.25);
[pks2, locs2] = findpeaks(bed2(end,xa:xb),dist2(xa:xb),'MinPeakDistance',0.25);
[dal2, locsd2] = findpeaks(-bed2(end,xa:xb),dist2(xa:xb),'MinPeakDistance',0.25);
[pks3, locs3] = findpeaks(bed3(end,xa:xb),dist3(xa:xb),'MinPeakDistance',0.25);
[dal3, locsd3] = findpeaks(-bed3(end,xa:xb),dist3(xa:xb),'MinPeakDistance',0.26);
[pks4, locs4] = findpeaks(bed4(end,xa:xb),dist4(xa:xb),'MinPeakDistance',0.25);
[dal4, locsd4] = findpeaks(-bed4(end,xa:xb),dist4(xa:xb),'MinPeakDistance',0.25);

[pks10, locs10] = findpeaks(-depth10(xa:xb,2)',dist1(xa:xb),'MinPeakDistance',0.28);
[dal10, locsd10] = findpeaks(depth10(xa:xb,2)',dist1(xa:xb),'MinPeakDistance',0.25);
[pks00, locs00] = findpeaks(-depth(xa:xb,2)',dist1(xa:xb),'MinPeakDistance',0.25);
[dal00, locsd00] = findpeaks(depth(xa:xb,2)',dist1(xa:xb),'MinPeakDistance',0.25);

% idx_pks=[1:10 12];
% idx_dal=[1:10];

locsd3=locsd3(1:end-1);
dal3=dal3(1:end-1);
locsd4=locsd4(1:end-1);
dal4=dal4(1:end-1);
locsd10=locsd10(1:end-1);
dal10=dal10(1:end-1);
locs10=locs10([1:10 12]);
pks10=pks10([1:10 12]);
%%
figure;
plot(dist1(xa:xb),-depth10(xa:xb,2))
hold on
scatter(locs10,pks10)
scatter(locsd10,-dal10)
hold on
plot(dist4(xa:xb),bed4(end,xa:xb))
hold on
scatter(locs4,pks4)
scatter(locsd4,-dal4)
%%
RMSEpiek1 = RootMeanSquareError(pks1,pks10,length(pks10));
RMSEdal1 = RootMeanSquareError(dal1,dal10,length(dal10));
RMSElocpeak1 = RootMeanSquareError(locs1,locs10,length(locs10));
RMSElocdal1 = RootMeanSquareError(locsd1,locsd10,length(locsd10));

RMSEpiek2 = RootMeanSquareError(pks2,pks10,length(pks10));
RMSEdal2 = RootMeanSquareError(dal2,dal10,length(dal10));
RMSElocpeak2 = RootMeanSquareError(locs2,locs10,length(locs10));
RMSElocdal2 = RootMeanSquareError(locsd2,locsd10,length(locsd10));

RMSEpiek3 = RootMeanSquareError(pks3,pks10,length(pks10));
RMSEdal3 = RootMeanSquareError(dal3,dal10,length(dal10));
RMSElocpeak3 = RootMeanSquareError(locs3,locs10,length(locs10));
RMSElocdal3 = RootMeanSquareError(locsd3,locsd10,length(locsd10));

RMSEpiek4 = RootMeanSquareError(pks4,pks10,length(pks10));
RMSEdal4 = RootMeanSquareError(dal4,dal10,length(dal10));
RMSElocpeak4 = RootMeanSquareError(locs4,locs10,length(locs10));
RMSElocdal4 = RootMeanSquareError(locsd4,locsd10,length(locsd10));

%%
mig_crest_data = nanmean((locs10-locs00)*100); %% x 1000 voor km naar m en gedeel door 10 . voor 10 jaar -> m/year
mig_trough_data = nanmean((locsd10-locsd00)*100); %% x 1000 voor km naar m en gedeel door 10 . voor 10 jaar -> m/year

mig_crest_RR = nanmean((locs1-locs00)*100); %% x 1000 voor km naar m en gedeel door 10 . voor 10 jaar -> m/year
mig_trough_RR = nanmean((locsd1-locsd00)*100); %% x 1000 voor km naar m en gedeel door 10 . voor 10 jaar -> m/year

mig_crest_RZ = nanmean((locs2-locs00)*100); %% x 1000 voor km naar m en gedeel door 10 . voor 10 jaar -> m/year
mig_trough_RZ = nanmean((locsd2-locsd00)*100); %% x 1000 voor km naar m en gedeel door 10 . voor 10 jaar -> m/year

mig_crest_UZ = nanmean((locs3-locs00)*100); %% x 1000 voor km naar m en gedeel door 10 . voor 10 jaar -> m/year
mig_trough_UZ = nanmean((locsd3-locsd00)*100); %% x 1000 voor km naar m en gedeel door 10 . voor 10 jaar -> m/year

mig_crest_ZZ = nanmean((locs4-locs00)*100); %% x 1000 voor km naar m en gedeel door 10 . voor 10 jaar -> m/year
mig_trough_ZZ = nanmean((locsd4-locsd00)*100); %% x 1000 voor km naar m en gedeel door 10 . voor 10 jaar -> m/year

mig_data = (mig_crest_data+mig_trough_data)/2;
mig_RR = (mig_crest_RR+mig_trough_RR)/2;
mig_RZ = (mig_crest_RZ+mig_trough_RZ)/2;
mig_UZ = (mig_crest_UZ+mig_trough_UZ)/2;
mig_ZZ = (mig_crest_ZZ+mig_trough_ZZ)/2;