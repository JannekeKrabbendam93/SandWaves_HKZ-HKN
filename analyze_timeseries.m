clear all
close all
clc

dt = 10; % min

addpath /Users/janneke/Documents/openearthtools2/openearthtools/trunk/matlab 
oetsettings

dir = '/Volumes/SEAGATE_EXP/2018-2019/Marid/';
file = 'RZall.bct';

workdir = ['/Volumes/SEAGATE_EXP/2018-2019/Marid/Componenten/Hydro_RZcomp/'];
%%
bct_basis = bct_io('read',[dir,file]);

Rz=bct_basis.Table(1).Data(:,2);
Zn=bct_basis.Table(2).Data(:,2);
time = bct_basis.Table(1).Data(:,1)*60-bct_basis.Table(1).Data(1,1)*60;

tstartz = [736960.534722222];
tstartn = [736795];
dtime = 0.1667*24;
tstart = datenum(2017,12,02,09,40,00);
tend=datenum(2018,02,02,12,40,00);
tall=tstart:10/1440:tend;

%%


%%
trz=t_tide(Rz,'interval',diff(time)/3600,'start time' ,tstart,'latitude',52.29472);

im2z=all(ismember(trz.name,'M2 '),2);
im4z=all(ismember(trz.name,'M4 '),2);
is2z=all(ismember(trz.name,'S2 '),2);
il2z=all(ismember(trz.name,'L2 '),2);
in2z=all(ismember(trz.name,'N2 '),2);
ims4z=all(ismember(trz.name,'MS4 '),2);
im6z=all(ismember(trz.name,'M6 '),2);
im8z=all(ismember(trz.name,'M8 '),2);
imu2z=all(ismember(trz.name,'MU2 '),2);
ik1z=all(ismember(trz.name,'K1 '),2);
io1z=all(ismember(trz.name,'O1 '),2);
iq1z=all(ismember(trz.name,'Q1 '),2);
is4z=all(ismember(trz.name,'S4 '),2);

im2z(1)=0;
im4z(1)=0;
im8z(1)=0;
im6z(1)=0;
im4sz(1)=0;
imu2z(1)=0;

idxz=logical(im2z+im4z+is2z+in2z+il2z+ims4z+im6z+im8z+imu2z+ik1z+io1z+iq1z+is4z);
trzM2.name=trz.name(idxz,:);
trzM2.freq=trz.freq(idxz);
trzM2.tidecon=trz.tidecon(idxz,:);
trzM2.type=trz.type;
trzM2.z0=trz.z0;
trzM2.dz0=trz.dz0;
trzM2.lat=trz.lat;
trzM2.period=trz.period;

%%
% [~,~,TIDECON,XOUT]=t_tide(Rz,'interval',ones(length(Rz)-1,1)*dtime,'start time' ,tstart,'latitude',52.29472,'sort','<->snr');
%%

%%
Rzpredic=t_predic(tall,trzM2,'latitude',52.29472,'synthesis',0);

%%
tzn=t_tide(Zn,'interval',diff(time)/3600,'start time' ,tstart,'latitude',52.29472);

im2n=all(ismember(tzn.name,'M2 '),2);
im4n=all(ismember(tzn.name,'M4 '),2);
is2n=all(ismember(tzn.name,'S2 '),2);
il2n=all(ismember(tzn.name,'L2 '),2);
in2n=all(ismember(tzn.name,'N2 '),2);
ims4n=all(ismember(tzn.name,'MS4 '),2);
im6n=all(ismember(tzn.name,'M6 '),2);
im8n=all(ismember(tzn.name,'M8 '),2);
imu2n=all(ismember(tzn.name,'MU2 '),2);
ik1n=all(ismember(tzn.name,'K1 '),2);
io1n=all(ismember(tzn.name,'O1 '),2);
iq1n=all(ismember(tzn.name,'Q1 '),2);
is4n=all(ismember(tzn.name,'S4 '),2);

im2n(1)=0;
im4n(1)=0;
im6n(1)=0;
im8n(1)=0;
ims4n(1)=0;
imu2n(1)=0;
idxn=logical(im2n+im4n+is2n+in2n+il2n+ims4n+im6n+im8n+imu2n+ik1n+io1n+iq1n+is4n);
tznM2.name=tzn.name(idxn,:);
tznM2.freq=tzn.freq(idxn);
tznM2.tidecon=tzn.tidecon(idxn,:);
tznM2.type=tzn.type;
tznM2.z0=tzn.z0;
tznM2.dz0=tzn.dz0;
tznM2.lat=tzn.lat;
tznM2.period=tzn.period;

%%
Znpredic=t_predic(tall,tznM2,'latitude',52.29472,'synthesis',0);

RMSE_Rz2 = RootMeanSquareError(Rzpredic,Rz(1:length(tall)),length(tall));
RMSE_Zn2 = RootMeanSquareError(Znpredic,Zn(1:length(tall)),length(tall));

%%
rest_rz = Rz(1:length(tall))-Rzpredic';
rest_zn = Zn(1:length(tall))-Znpredic';

%%
figure;
plot(time(1:length(tall))/3600/24,rest_rz,'LineWidth',1.5)
set(gca,'FontSize',22,'FontName','Times New Roman')
xlabel('t [days]')
title('Rest Rz')
xlim([0 62])
ylim([-1 1.5])
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1400,400]);
set(gcf,'PaperPositionMode','auto')
%%
figure;
plot(time(1:length(tall))/3600/24,rest_zn,'LineWidth',1.5)
set(gca,'FontSize',22,'FontName','Times New Roman')
xlabel('t [days]')
title('Rest Zn')
xlim([0 62])
ylim([-1 1.5])
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1400,400]);
set(gcf,'PaperPositionMode','auto')
%%
figure;
plot(time/3600/24,Rz,'LineWidth',1.5)
hold on
plot(time(1:length(tall))/3600/24,Rzpredic,'LineWidth',1.5)
set(gca,'FontSize',22,'FontName','Times New Roman')
xlabel('t [days]')
ylabel('R (t) [m/s]')
h=legend('data','sum 13 components')
set(h.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));
xlim([0 62])
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1400,400]);
set(gcf,'PaperPositionMode','auto')

%%
figure;
plot(time/3600/24,Zn,'LineWidth',1.5)
hold on
% plot(time/3600,Zncomp)
plot(time(1:length(tall))/3600/24,Znpredic,'LineWidth',1.5)
set(gca,'FontSize',22,'FontName','Times New Roman')
xlabel('t [days]')
ylabel('Z (t) [m]')
h=legend('data','sum 13 components')
set(h.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.7]));
xlim([0 62])
ylim([-1.5 2])
set(gca,'GridAlpha',0.5);
grid on; 
grid minor
box on
set(gcf,'Position',[100, 100,1400,400]);
set(gcf,'PaperPositionMode','auto')

%%
trz_rest=t_tide(rest_rz,'interval',diff(time(1:length(rest_zn)))/3600,'start time' ,tstart,'latitude',52.29472,'sort','<->snr');
% %%
% w(1) = 2.*pi/((360./28.984104).*3600);                                        % M2
% w(2) = 2.*pi/((360./30).*3600);                                               % S2
% w(3) = 2.*pi/((360./28.43973).*3600);                                         % N2
% w(4) = 2.*pi/((360./57.96821).*3600);                                         % M4
% w(5) = 2.*pi/((360./86.95232).*3600);                                         % M6
% w(6) = 2.*pi/((360./115.93642).*3600);                                       % M8
% w(7) = 2.*pi/((360./15.041069).*3600);                                        % K1
% w(8) = 2.*pi/((360./13.943035).*3600);                                        % O1
% w(9) = 2.*pi/((360./44.025173).*3600);                                        % MK3
% w(10) = 2.*pi/((360./60).*3600);                                               % S4 (not prescribed in kuststrook but still included in analysis)
% w(11) = 2.*pi/((360./57.423832).*3600);                                       % MN4
% w(12) = 2.*pi/((360./58.984104).*3600);                                       % MS4
% w(13) = 2.*pi/((360./28.512583).*3600);                                       % NU2
% w(14) = 2.*pi/((360./90).*3600);                                              % S6 (not prescribed in kustrook but still included in analysis)
% w(15) = 2.*pi/((360./27.968208).*3600);                                       % MU2
% w(16) = 2.*pi/((360./13.398661).*3600);                                       % Q1
% w(17) = 2.*pi/((360./29.958933).*3600);                                       % T2
% %R2 (not prescribed in kuststrook)
% %2Q1
% w(18) = 2.*pi/((360./14.958931).*3600);                                       % P1
% w(19) = 2.*pi/((360./31.015896).*3600);                                       % 2SM2
% %M3 (not prescribed in kuststrook)
% w(20) = 2.*pi/((360./29.5284789).*3600);                                        % L2
% w(21) = 2.*pi/((360./42.92714).*3600);                                        % 2MK3
% w(22) = 2.*pi/((360./30.082138).*3600);                                        % K2
% 
% % %%
% % dph=45;
% % 
% % testfreq=trz.freq.*2*pi/3600;
% % testamp=trz.tidecon(:,1);
% % testphase=(dph+trz.tidecon(:,3)).*pi/180;
% % 
% % for aa=1:length(time)
% %     RZallt(aa)=trz.z0+sum(testamp.*cos(testfreq.*time(aa)-testphase));
% % end
% %%
% dph=45;
% Rzcomp=0.0556+0.754*cos(w(1)*time-(dph+75.55)*pi/180)+0.149*cos(w(4)*time-(dph+105.27)*pi/180)...   % M0+M2+M4
%     +0.195*cos(w(2)*time-(dph+127.41)*pi/180)+0.09*cos(w(12)*time-(dph+152.36)*pi/180)...           % S2+MS4
%     +0.139*cos(w(3)*time-(55.21)*pi/180)+0.066*cos(w(20)*time-(86.03)*pi/180);%...             % N2+L2
%    % +0.0607*cos(w(7)*time-9.20*pi/180)+0.0709*cos(w(8)*time-182.92*pi/180)...             % K1+O1
%   %  +0.0165*cos(w(16)*time-137.45*pi/180)+0.056*cos(w(11)*time-100.42*pi/180)...       %Q1+MN4
%  %  +0.0522*cos(w(15)*time-185.79*pi/180)+0.03*cos(w(5)*time-132.02*pi/180)...%MU2+M6
%  %  +0.0131*cos(w(6)*time-189.92*pi/180)+0.0064 *cos(w(10)*time-238.64*pi/180);       %M8+S4
% dph2=45;
% Zncomp=0.0002+0.4862*cos(w(1)*time-(dph2+127.13)*pi/180)+0.1582*cos(w(4)*time-(142.30+dph2)*pi/180)...    % M0+M2+M4
%     +0.1224*cos(w(2)*time-(dph2+201.19)*pi/180)+0.0738*cos(w(12)*time-(dph2+204.48)*pi/180)...            % S2+MS4
%     +0.0757*cos(w(3)*time-(dph2+110.17)*pi/180)+0.0524*cos(w(20)*time-(dph2+128.51)*pi/180);%...             % N2+L2
%     %+0.0953*cos(w(7)*time-342.81*pi/180)+0.113*cos(w(8)*time-176.51*pi/180)...                 % K1+O1
%     %+0.0432*cos(w(16)*time-136.53*pi/180)+0.0504*cos(w(11)*time-113.88*pi/180)...       %Q1+MN4
%     %+0.066*cos(w(15)*time-213.57*pi/180)+0.058*cos(w(5)*time-256.38*pi/180);       %MU2+M6
% %%
% RMSE_Rz = RootMeanSquareError(Rzcomp,Rz,length(time));
% RMSE_Zn = RootMeanSquareError(Zncomp,Zn,length(time));
% 
% % %%
% % trimfile1=fullfile(workdir,'trim-bas_sed.dat');
% % trim1 = vs_use(trimfile1);
% % m1 = qpfopen(trim1.DatExt);
% % l1 = qpread(m1);
% % 
% % tt=vs_time(trim1);
% % Morfotime=tt.t/60;
% % ttime=Morfotime*60;
% % 
% % U1 = vs_get(trim1,'map-series',{1:length(Morfotime)},'U1',{2 2 0});
% % wl1 = vs_get(trim1,'map-series',{1:length(Morfotime)},'S1',{2 2});
% 
% % %%
% % Uz = nanmean(cell2mat(U1),3);
% % Zz = squeeze(cell2mat(wl1));
% % 
% % Rzmod = Uz+sqrt(9.81/22.5)*Zz;