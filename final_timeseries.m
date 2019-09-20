close all
clear all
clc

dt = 10; % min

addpath /Users/janneke/Documents/openearthtools2/openearthtools/trunk/matlab 
oetsettings

%% HKN data - Water level
Z1N = xlsread('/Users/janneke/Documents/2018-2019/Data/WL_HKN/December_2017/02-Final/01 - data/HKNB/HKN_20180409_Fugro_MetOcean Buoys HKNB December 2017 WaterLevelNonInterpolated_F.xlsx');
Z2N = xlsread('/Users/janneke/Documents/2018-2019/Data/WL_HKN/January_2018/03-Final/01 - Data/HKNB/HKN_20180412_Fugro_MetOcean Buoys HKNB January 2018 WaterLevelNonInterpolated_F.xlsx');
Z3N = xlsread('/Users/janneke/Documents/2018-2019/Data/WL_HKN/February_2018/02-Final/01 - Data/HKNB/HKN_20180416_Fugro_MetOcean Buoys HKNB February 2018 WaterLevelNonInterpolated_F.xlsx');

%% HKN data - Current
C1N = xlsread('/Users/janneke/Documents/2018-2019/Data/WL_HKN/December_2017/02-Final/01 - data/HKNB/HKN_20180409_Fugro_MetOcean Buoys HKNB December 2017 WaveCurrentDataStat_F.xlsx');
C2N = xlsread('/Users/janneke/Documents/2018-2019/Data/WL_HKN/January_2018/03-Final/01 - Data/HKNB/HKN_20180412_Fugro_MetOcean Buoys HKNB January 2018 WaveCurrentDataStat_F.xlsx');
C3N = xlsread('/Users/janneke/Documents/2018-2019/Data/WL_HKN/February_2018/02-Final/01 - Data/HKNB/HKN_20180416_Fugro_MetOcean Buoys HKNB February 2018 WaveCurrentDataStat_F.xlsx');

%% HKZ data - Water level
Z1Z = xlsread('/Users/janneke/Documents/2018-2019/Data/WL_HKZ/December_2017/02-Final/01 - Data/HKZA/HKZ_20180207_Fugro_MetOcean Buoys HKZA December 2017 WaterLevelNonInterpolated_F.xlsx');
Z2Z = xlsread('/Users/janneke/Documents/2018-2019/Data/WL_HKZ/January_2018/02-Final/01 - data/HKZA/HKZ_20180412_Fugro_MetOcean Buoys HKZA January 2018 WaterLevelNonInterpolated_F.xlsx');
Z3Z = xlsread('/Users/janneke/Documents/2018-2019/Data/WL_HKZ/February_2018/02-Final/01 - Data/HKZA/HKZ_20180413_Fugro_MetOcean Buoys HKZA February 2018 WaterLevelNonInterpolated_F.xlsx');

%% HKZ data - Current
C1Z = xlsread('/Users/janneke/Documents/2018-2019/Data/WL_HKZ/December_2017/02-Final/01 - Data/HKZA/HKZ_20180207_Fugro_MetOcean Buoys HKZA December 2017 WaveCurrentDataStat_F.xlsx');
C2Z = xlsread('/Users/janneke/Documents/2018-2019/Data/WL_HKZ/January_2018/02-Final/01 - data/HKZA/HKZ_20180412_Fugro_MetOcean Buoys HKZA January 2018 WaveCurrentDataStat_F.xlsx');
C3Z = xlsread('/Users/janneke/Documents/2018-2019/Data/WL_HKZ/February_2018/02-Final/01 - Data/HKZA/HKZ_20180413_Fugro_MetOcean Buoys HKZA February 2018 WaveCurrentDataStat_F.xlsx');

%%
ZN_all = [Z1N(:,2); Z2N(:,2); Z3N(:,2)];
meanZN = nanmean(ZN_all);
ZN_nan = isnan(ZN_all);
tN = (1:numel(ZN_all))*dt*60;
tN2 = datenum(tN);

ZN_all(ZN_nan)=interp1(tN(~ZN_nan),ZN_all(~ZN_nan),tN(ZN_nan)); 
depN = ZN_all;
ZN_all = ZN_all - meanZN;

%%
dstep2 = (50*60)/10;   % 50 uur averaging interval, waarbij 10 de interval van de metingen is
F2 = fspecial('average',[dstep2 1]);    % averaging filter

nandepN = nanconv(depN, F2, 'edge','nanout');  % bepalen residual
depN2 = depN-nandepN;  % aftrekken van tijdserie
nanZN = nanconv(ZN_all, F2, 'edge','nanout');  % bepalen residual
ZN_all2 = ZN_all-nanZN;  % aftrekken van tijdserie
M0= 0.05;
%% Interpolate water level data - HKZ
ZZ_all = [Z1Z(:,2); Z2Z(:,2);Z3Z(:,2)]/1E6;
meanZZ = nanmean(ZZ_all);
ZZ_nan = isnan(ZZ_all);
tZ = (1:numel(ZZ_all))*dt*60;

ZZ_all(ZZ_nan)=interp1(tZ(~ZZ_nan),ZZ_all(~ZZ_nan),tZ(ZZ_nan)); 
depZ = ZZ_all;
ZZ_all = ZZ_all - meanZZ;

nandepZ = nanconv(depZ, F2, 'edge','nanout');  % bepalen residual
depZ2 = depZ-nandepZ;  % aftrekken van tijdserie
nanZZ = nanconv(ZZ_all, F2, 'edge','nanout');  % bepalen residual
ZZ_all2 = ZZ_all-nanZZ;  % aftrekken van tijdserie

%% Create time series of direction (T) and magnitude (M) of current at all depths - HKN
depth = 4:2:30;
select = 4:17;
for i = 1:length(select)    
    j = select(i);
    DN_all(i,:) = [C1N(:,j); C2N(:,j); C3N(:,j)];
end

select2 = 18:31;
for m=1:length(select2)
    o = select2(m);
    MN_all(m,:) = [C1N(:,o); C2N(:,o); C3N(:,o)]./100;
end

%% Create time series of direction (T) and magnitude (M) of current at all depths - HKZ
for i = 1:length(select)    
    j = select(i);
    DZ_all(i,:) = [C1Z(:,j); C2Z(:,j); C3Z(:,j)]./1E6;
end

for m=1:length(select2)
    o = select2(m);
    MZ_all(m,:) = [C1Z(:,o); C2Z(:,o); C3Z(:,o)]./1E8;
end

%% Calculate U and V from magnitude and direction in direction of grid - HKN
for i = 1:length(depth)
    for m = 1:length(DN_all(i,:))
        UN1(i,m) = sind(DN_all(i,m))*MN_all(i,m);
        VN1(i,m) = cosd(DN_all(i,m))*MN_all(i,m);
    end
end

for i = 1:length(depth)
    for c = 1:length(DZ_all(i,:))
        UZ1(i,c) = sind(DZ_all(i,c))*MZ_all(i,c);
        VZ1(i,c) = cosd(DZ_all(i,c))*MZ_all(i,c);
    end
end

%%
UaveZ=nanmean(UZ1(1:9,:),1);
VaveZ=nanmean(VZ1(1:9,:),1);
UaveN=nanmean(UN1(1:9,:),1);
VaveN=nanmean(VN1(1:9,:),1);

%%
figure;
scatter(UaveZ,VaveZ,'filled','MarkerEdgeColor',[0 0 0])
set(gca,'FontSize',22,'FontName','Times New Roman')
xlabel('U (m s^{-1})')
ylabel('V (m s^{-1})')
axis equal
ylim([-1 1])
xlim([-1 1])
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,600,600]);
set(gcf,'PaperPositionMode','auto')
eval(['print -depsc -tiff ' ,['UV_hkz.eps']])
%%
figure;
scatter(UaveN,VaveN,'filled','MarkerEdgeColor',[0 0 0])
set(gca,'FontSize',22,'FontName','Times New Roman')
xlabel('U (m s^{-1})')
ylabel('V (m s^{-1})')
axis equal
ylim([-1 1])
xlim([-1 1])
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,600,600]);
set(gcf,'PaperPositionMode','auto')
eval(['print -depsc -tiff ' ,['UV_hkn.eps']])
%%
UZ2 = UaveZ*sind(71)+VaveZ*cosd(71);
UN2 = UaveN*sind(71)+VaveN*cosd(71);

VZ2 = -UaveZ*cosd(71)+VaveZ*sind(71);
VN2 = -UaveN*cosd(71)+VaveN*sind(71);

figure;
scatter(UaveZ,VaveZ,'filled','MarkerEdgeColor',[0 0 0])
hold on
scatter(UZ2,VZ2,'filled','MarkerEdgeColor',[0 0 0])
set(gca,'FontSize',22,'FontName','Times New Roman')
xlabel('U (m s^{-1})')
ylabel('V (m s^{-1})')
axis equal
ylim([-1 1])
xlim([-1 1])
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,600,600]);
set(gcf,'PaperPositionMode','auto')
eval(['print -depsc -tiff ' ,['UV_hkzrot.eps']])
%%
figure;
scatter(UaveN,VaveN,'filled','MarkerEdgeColor',[0 0 0])
hold on
scatter(UN2,VN2,'filled','MarkerEdgeColor',[0 0 0])
set(gca,'FontSize',22,'FontName','Times New Roman')
xlabel('U (m s^{-1})')
ylabel('V (m s^{-1})')
axis equal
ylim([-1 1])
xlim([-1 1])
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,600,600]);
set(gcf,'PaperPositionMode','auto')
eval(['print -depsc -tiff ' ,['UV_hknrot.eps']])
%%
figure;
plot(UN2)
hold on
plot(UZ2)

%%
figure;
plot(VN2)
hold on
plot(VZ2)
%%
UZ_nan = isnan(UZ2);
tuZ = (1:numel(UZ2))*dt*60;
UZ2(UZ_nan)=interp1(tuZ(~UZ_nan),UZ2(~UZ_nan),tuZ(UZ_nan)); 

UN_nan = isnan(UN2);
tuN = (1:numel(UN2))*dt*60;
UN2(UN_nan)=interp1(tuN(~UN_nan),UN2(~UN_nan),tuN(UN_nan)); 

%%
dstep3 = (50*60)/10;   % 75 uur averaging interval, waarbij 10 de interval van de metingen is
F3 = fspecial('average',[dstep3]);    % averaging filter
%%
nanuN = nanconv(UN2, F3, 'edge','nanout');  % bepalen residual
UN2 = UN2-0.9*nanuN;  % aftrekken van tijdserie
nanuZ = nanconv(UZ2, F3, 'edge','nanout');  % bepalen residual
UZ2 = UZ2-0.3*nanuZ;  % aftrekken van tijdserie
nanvN = nanconv(VN2, F3, 'edge','nanout');  % bepalen residual
VN2 = VN2-nanvN;  % aftrekken van tijdserie
nanvZ = nanconv(VZ2, F3, 'edge','nanout');  % bepalen residual
VZ2 = VZ2-nanvZ;  % aftrekken van tijdserie

%%
figure;
plot(UZ2)
hold on
plot(UZ2)
%%
Daystart = datetime(2017,12,02,09,40,00);
Dayend = datetime(2018,02,02,20,10,00);
Timedays = Daystart:(10/1440):Dayend;
length(Timedays)
time = datenum(Timedays);

DaystartU = datetime(2017,12,01,00,00,00);
DayendU = datetime(2018,02,28,23,50,00);
TimedaysU = DaystartU:(10/1440):DayendU;
length(TimedaysU)
timeU = datenum(TimedaysU);

DaystartZ = datetime(2017,12,01,00,00,00);
DayendZ = datetime(2018,02,02,20,10,00);
TimedaysZ = DaystartZ:(10/1440):DayendZ;
length(TimedaysZ)
timeZ = datenum(TimedaysZ);

DaystartN = datetime(2017,12,02,09,40,00);
DayendN = datetime(2018,02,28,23,50,00);
TimedaysN = DaystartN:(10/1440):DayendN;
length(TimedaysN)
timeN = datenum(TimedaysN);

%%
t1u = find(TimedaysU==Timedays(1));
t2u = find(TimedaysU==Timedays(end));

t1n = find(TimedaysN==Timedays(1));
t2n = find(TimedaysN==Timedays(end));

t1z = find(TimedaysZ==Timedays(1));
t2z = find(TimedaysZ==Timedays(end));

%%
depZ3=depZ2+meanZZ;
depN3=depN2+meanZN;

%%
for aa=1:length(depZ3);
    tempZ(aa)=sqrt(9.81/depZ3(aa))*ZZ_all2(aa);
end

for bb=1:length(depN3);
    tempN(bb)=sqrt(9.81/depN3(bb)).*ZN_all2(bb);
end
%%
Rz=smooth(UZ2(t1u:t2u)+tempZ(t1z:t2z));
Rn=smooth(UN2(t1u:t2u)-tempN(t1n:t2n));
Uz=smooth(UZ2(t1u:t2u));
Un=smooth(UN2(t1u:t2u));
Zz=smooth(ZZ_all2(t1z:t2z)');
Zn=smooth(ZN_all2(t1n:t2n)');

save('buoy_data.mat','Uz','Un','Zn','Zz','UaveZ','VaveZ','UaveN','VaveN')
%%
dir = '/Volumes/SEAGATE_EXP/2018-2019/Marid/Timeseries/RRtime_hc_MF50/';
file = 'RR_MF50.bct';

%%
bct_basis = bct_io('read',[dir,file]);
Tstart = 2020;
T100 = 5.2603E4;
T50 = 1.0521E5;
T60 = 8.7690E4;
T25 = 2.1044E5;
Tspin = 3000;
T50_40yr = 4.2087E5;

Tall = length(time)*10+Tstart-1;
Tstop = Tstart+Tspin+T25;
Time=Tstart:10:Tstop;

%%
Rz2=[Rz(1:end-40); Rz(1:end-40); Rz(1:end-40); Rz(1:end-40); Rz];
Uz2=[Uz(1:end-43); Uz(1:end-43); Uz(1:end-43); Uz(1:end-43); Uz];
Un2=[Un(1:end-43); Un(1:end-43); Un(1:end-43); Un(1:end-43); Un];
Rn2=[Rn(1:end-43); Rn(1:end-43); Rn(1:end-43); Rn(1:end-43); Rn];
Zz2=[Zz(1:end-60); Zz(1:end-60); Zz(1:end-60); Zz(1:end-60); Zz];
Zn2=[Zn(1:end-60); Zn(1:end-60); Zn(1:end-60); Zn(1:end-60); Zn];

%%
figure;
plot(Rz(length(Zn)/2-500:length(Rz)/2+500))

%%
bct_final = bct_basis;

lala=[Time' Rz2(1:length(Time)) Rz2(1:length(Time))];
lolo=[Time' Rn2(1:length(Time)) Rn2(1:length(Time))];
bct_final.Table(1).Data=lala;
bct_final.Table(2).Data=lolo;

Rzfin=Rz2(1:length(Time));
Rnfin=Rn2(1:length(Time));
Uzfin=Uz2(1:length(Time));
Unfin=Un2(1:length(Time));
Znfin=Zn2(1:length(Time));
Zzfin=Zz2(1:length(Time));

figure;
plot(Rzfin)

nanmean(Rzfin)
nanmean(Rnfin)
nanmean(Znfin)
nanmean(Uz2)
nanmean(Un2)

save('bc_all.mat','Uzfin','Unfin','Zzfin','Znfin','Time')
%% 

bct_io('write','RR_MF25_2.bct',bct_final);
%%
dt = Time/24/60-Time(1)/24/60-2.08;

%%
figure;
plot(Rzfin)

figure;
plot(Rnfin)
%%
fig = figure;
left_color = [0 0 1];
right_color = [1 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(dt,Uzfin,'-b','LineWidth',2)
xlim([33 40])
set(gca,'FontSize',22,'FontName','Times New Roman')
ylabel('U [m/s]')
ylim([-1 1])
yyaxis right
plot(dt,Zzfin,'-r','LineWidth',2)
xlim([26 40])
ylim([-1.5 1.5])
set(gca,'FontSize',22,'FontName','Times New Roman')
title('HKZ')
xlabel('t [days]')
ylabel('\zeta [m]')


set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1400,400]);
set(gcf,'PaperPositionMode','auto')
eval(['print -depsc -tiff ' ,['UZ_hkz_day30-40.eps']])
%%
fig = figure;
left_color = [0 0 1];
right_color = [1 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(dt,Unfin,'-b','LineWidth',2)
xlim([26 40])
ylim([-1 1])
set(gca,'FontSize',22,'FontName','Times New Roman')
ylabel('U [m/s]')

yyaxis right
plot(dt,Znfin,'-r','LineWidth',2)
xlim([26 40])
ylim([-1.5 1.5])
set(gca,'FontSize',22,'FontName','Times New Roman')
title('HKN')
xlabel('t [days]')
ylabel('\zeta [m]')
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1400,400]);
set(gcf,'PaperPositionMode','auto')
eval(['print -depsc -tiff ' ,['UZ_hkn_day30-40.eps']])
%%
figure;
plot(Time/60/24-Time(1)/60/24,Un,'LineWidth',2)
set(gca,'FontSize',22,'FontName','Times New Roman')
title('HKN')
xlabel('t [days]')
ylabel('\zeta [m]')
% ylim([-4E-10 14E-10])
xlim([0 60])
% legend(['L=',num2str(L(1))],['L=',num2str(L(2))],['L=',num2str(L(3))],['L=',num2str(L(4))],['L=',num2str(L(5))],['L=',num2str(L(6))]);%,['L=',num2str(L(7))])
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1000,600]);
set(gcf,'PaperPositionMode','auto')
eval(['print -depsc -tiff ' ,['Zn.eps']])
% hold on
% plot(lolo(:,2))

