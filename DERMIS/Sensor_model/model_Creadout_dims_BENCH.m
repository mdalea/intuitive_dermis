% remember this only works for AC VTX, results are AC magnitudes
%Cs1 and Cs2 are reversed but who cares

%plots the envelope to the square/sine wave Vtx

clear all
close all

%% ==== SHEAR FORCE SWEEP, NO NORMAL FORCE =====
% C0=0.6125e-15; % no normal force
N=5;
% deltaC=[0:C0/2^N:C0];
Cf=40e-15; %10e-15;
% Cs1=C0-deltaC; Cs2=C0+deltaC;
VDD=3; %VDD=3.3;
Vocm=VDD/2; % this is the source of nonlinearity for VCM measurements
Coffset1=1.7e-12%0.85e-12 %1e-15;
Coffset2=1.7e-12%0.85e-12 %1e-15;

x0=50e-6; % overlap at force (shear) = 0
z0=125e-6; % thickness at force (normal) = 0
y0=125e-6; % constant height overlap
deltaX_max=50e-6;
deltaX_res=deltaX_max./2^N;
deltaX=[-deltaX_max:deltaX_res:deltaX_max];
deltaZ=0;

eps=8.85e-12;
eps_pdms=2.77*600; %factor 300x larger capacitance

Cs1a=(eps_pdms.*eps.*(x0-deltaX).*y0)./(z0-deltaZ);
Cs1b=(eps_pdms.*eps.*(x0-deltaX).*y0)./(z0-deltaZ);
Cs1=Coffset1+(Cs1a.*Cs1b)./(Cs1a+Cs1b);
%Cs1(isnan(Cs1))=0

Cs2a=(eps_pdms.*eps.*(x0+deltaX).*y0)./(z0-deltaZ);
Cs2b=(eps_pdms.*eps.*(x0+deltaX).*y0)./(z0-deltaZ);
Cs2=Coffset2+(Cs2a.*Cs2b)./(Cs2a+Cs2b);
%Cs2(isnan(Cs2))=0

deltaC=(eps_pdms.*eps.*(deltaX).*y0)./(z0-deltaZ)./2;
figure
plot(deltaX,Cs1); hold on;
plot(deltaX,Cs2)
title('FSHEAR=sweep; FNORM=0')
xlabel('Delta-X (m)')
ylabel('Capacitance (F)')
legend('Cs1','Cs2','Location','best')

% VDIFF mode
Vinp=VDD; %0 IF CT; VDD IF DT
Vinn=VDD; %0 IF CT; VDD IF DT

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;
%Vod can be approximated as 2*(Vinp-Vinn)*(Cs/Cf) if Cs=Cs1=Cs2
%Voutp can be approx. as (Vinp-Vinn)*(Cs+delC)/Cf
%Voutn can be approx. as -(Vinp-Vinn)*(Cs-delC)/Cf

figure

subplot(1,2,1)
plot(deltaX,Vod)
xlabel('Delta-X (m)')
ylabel('Voutput differential (Voutp-Voutn)')
title('FSHEAR=sweep; FNORM=0')

subplot(1,2,2)
plot(deltaX,Voutp); hold on
plot(deltaX,Voutn)
xlabel('Delta-X (m)')
ylabel('V')
legend('Voutp','Voutn','Location','best')
title('VDIFF MODE (VINP=VINN=VDD)')

Voutp(2)
Voutn(2)
Vod(2)
Vod(end)

% VCM mode
Vinp=VDD %0 %VDD;
Vinn=0 %VDD %0;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

figure

subplot(1,2,1)
plot(deltaX,Vod)
xlabel('Delta-X (m)')
ylabel('Voutput differential (Voutp-Voutn)')
title('FSHEAR=sweep; FNORM=0')

subplot(1,2,2)
plot(deltaX,Voutp); hold on
plot(deltaX,Voutn)
xlabel('Delta-X (m)')
ylabel('V')
legend('Voutp','Voutn','Location','best')
title('VCM MODE (VINP=VDD; VINN=0)')

% VDIFF/VCM mode
Vinp=VDD;
Vinn=VDD;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

Vinp=VDD;
Vinn=0;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vcm=Voutp-Voutn;


figure
subplot(1,2,1)
plot(deltaX,Vod./Vcm)
xlabel('Delta-X (m)')
ylabel('Voutput ratio (Vod/Vocm)')
title('FSHEAR=sweep; FNORM=0')

subplot(1,2,2)
plot(deltaX,Vod); hold on
plot(deltaX,Vcm)
xlabel('Delta-X (m)')
ylabel('V')
legend('Vod','Vocm','Location','best')
title('VDIFF/VCM RATIOMETRIC MODE')
%% ==== SHEAR FORCE SWEEP; MID NORMAL FORCE =====

x0=50e-6; % overlap at force (shear) = 0
z0=125e-6 %for deltaZ_max=100e-6   %125e-6; % thickness at force (normal) = 0
y0=125e-6; % constant height overlap
deltaX_max=50e-6;
deltaX_res=deltaX_max./2^N;
deltaX=[-deltaX_max:deltaX_res:deltaX_max];
deltaZ=50e-6;

eps=8.85e-12;
eps_pdms=2.77;
Cs1a=(eps_pdms.*eps.*(x0-deltaX).*y0)./(z0-deltaZ);
Cs1b=(eps_pdms.*eps.*(x0-deltaX).*y0)./(z0-deltaZ);
Cs1=(Cs1a.*Cs1b)./(Cs1a+Cs1b);
Cs1(isnan(Cs1))=0

Cs2a=(eps_pdms.*eps.*(x0+deltaX).*y0)./(z0-deltaZ);
Cs2b=(eps_pdms.*eps.*(x0+deltaX).*y0)./(z0-deltaZ);
Cs2=(Cs2a.*Cs2b)./(Cs2a+Cs2b);

deltaC=(eps_pdms.*eps.*(deltaX).*y0)./(z0-deltaZ)./2;
figure
plot(deltaX,Cs1); hold on;
plot(deltaX,Cs2)
title('FSHEAR=sweep; FNORM=MID')
xlabel('Delta-X (m)')
ylabel('Capacitance (F)')
legend('Cs1','Cs2','Location','best')

% VDIFF mode
Vinp=VDD; %0 IF CT; VDD IF DT
Vinn=VDD; %0 IF CT; VDD IF DT

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

figure

subplot(1,2,1)
plot(deltaX,Vod)
xlabel('Delta-X (m)')
ylabel('Voutput differential (Voutp-Voutn)')
title('FSHEAR=sweep; FNORM=MID')

subplot(1,2,2)
plot(deltaX,Voutp); hold on
plot(deltaX,Voutn)
xlabel('Delta-X (m)')
ylabel('V')
legend('Voutp','Voutn','Location','best')
title('VDIFF MODE (VINP=VINN=VDD)')

Voutp(2)
Voutn(2)
Vod(2)
Vod(end)

% VCM mode
Vinp=VDD %0 %VDD;
Vinn=0 %VDD %0;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

figure

subplot(1,2,1)
plot(deltaX,Vod)
xlabel('Delta-X (m)')
ylabel('Voutput differential (Voutp-Voutn)')
title('FSHEAR=sweep; FNORM=MID')

subplot(1,2,2)
plot(deltaX,Voutp); hold on
plot(deltaX,Voutn)
xlabel('Delta-X (m)')
ylabel('V')
legend('Voutp','Voutn','Location','best')
title('VCM MODE (VINP=VDD; VINN=0)')

% VDIFF/VCM mode
Vinp=VDD;
Vinn=VDD;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

Vinp=VDD;
Vinn=0;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vcm=Voutp-Voutn;


figure
subplot(1,2,1)
plot(deltaX,Vod./Vcm)
xlabel('Delta-X (m)')
ylabel('Voutput ratio (Vod/Vocm)')
title('FSHEAR=sweep; FNORM=MID')

subplot(1,2,2)
plot(deltaX,Vod); hold on
plot(deltaX,Vcm)
xlabel('Delta-X (m)')
ylabel('V')
legend('Vod','Vocm','Location','best')
title('VDIFF/VCM RATIOMETRIC MODE')




%% === NO SHEAR FORCE; NORMAL FORCE ONLY

x0=50e-6; % overlap at force (shear) = 0
z0=125e-6; % thickness at force (normal) = 0
y0=125e-6; % constant height overlap
deltaZ_max=100e-6; %deltaZ_max=50e-6;
deltaZ_res=deltaZ_max./2^N;
deltaZ=[0:deltaZ_res:deltaZ_max];
%deltaZ=[-deltaZ_max:deltaZ_res:deltaZ_max];
deltaX=0;

eps=8.85e-12;
eps_pdms=2.77;
Cs1a=(eps_pdms.*eps.*(x0-deltaX).*y0)./(z0-deltaZ);
Cs1b=(eps_pdms.*eps.*(x0-deltaX).*y0)./(z0-deltaZ);
Cs1=(Cs1a.*Cs1b)./(Cs1a+Cs1b);
Cs1(isnan(Cs1))=0

Cs2a=(eps_pdms.*eps.*(x0+deltaX).*y0)./(z0-deltaZ);
Cs2b=(eps_pdms.*eps.*(x0+deltaX).*y0)./(z0-deltaZ);
Cs2=(Cs2a.*Cs2b)./(Cs2a+Cs2b);

deltaC=(eps_pdms.*eps.*(x0).*y0)./(deltaZ)./2;
figure
plot(deltaZ,Cs1); hold on;
plot(deltaZ,Cs2)
xlabel('Delta-Z (m)')
ylabel('Capacitance (F)')
legend('Cs1','Cs2','Location','best')
title('FSHEAR=0; FNORM=sweep')

% VDIFF mode
Vinp=VDD;
Vinn=VDD;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

figure
subplot(1,2,1)
plot(deltaZ,Vod)
xlabel('Delta-Z (m)')
ylabel('Voutput differential (Voutp-Voutn)')
title('FSHEAR=0; FNORM=sweep')

subplot(1,2,2)
plot(deltaZ,Voutp); hold on
plot(deltaZ,Voutn)
xlabel('Delta-Z (m)')
ylabel('V')
legend('Voutp','Voutn','Location','best')
title('VDIFF MODE (VINP=VINN=VDD)')

Voutp(2)
Voutn(2)
Vod(2)
Vod(end)

% VCM mode
Vinp=VDD;
Vinn=0;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

disp(['MAX Vod during VCM mode: ',num2str(max(Vod))])

figure
subplot(1,2,1)
plot(deltaZ,Vod)
xlabel('Delta-Z (m)')
ylabel('Voutput differential (Voutp-Voutn)')
title('FSHEAR=0; FNORM=sweep')

subplot(1,2,2)
plot(deltaZ,Voutp); hold on
plot(deltaZ,Voutn)
xlabel('Delta-Z (m)')
ylabel('V')
legend('Voutp','Voutn','Location','best')
title('VCM MODE (VINP=VDD; VINN=0)')

% VDIFF/VCM mode
Vinp=VDD;
Vinn=VDD;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

Vinp=VDD;
Vinn=0;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vcm=Voutp-Voutn;

Vratio=Vod./Vcm;

figure
subplot(1,2,1)
plot(deltaZ,Vratio)
xlabel('Delta-Z (m)')
ylabel('Voutput ratio (Vod/Vocm)')
title('FSHEAR=0; FNORM=sweep')

subplot(1,2,2)
plot(deltaZ,Vod); hold on
plot(deltaZ,Vcm)
xlabel('Delta-Z (m)')
ylabel('V')
legend('Vod','Vocm','Location','best')
title('VDIFF/VCM RATIOMETRIC MODE')

%% === MAX CONSTANT SHEAR FORCE; NORMAL FORCE ONLY

x0=50e-6; % overlap at force (shear) = 0
z0=125e-6; % thickness at force (normal) = 0
y0=125e-6; % constant height overlap
deltaZ_max=100e-6; %deltaZ_max=50e-6;
deltaZ_res=deltaZ_max./2^N;
deltaZ=[0:deltaZ_res:deltaZ_max];
deltaX=50e-6; %max

eps=8.85e-12;
eps_pdms=2.77;
Cs1a=(eps_pdms.*eps.*(x0-deltaX).*y0)./(z0-deltaZ);
Cs1b=(eps_pdms.*eps.*(x0-deltaX).*y0)./(z0-deltaZ);
Cs1=(Cs1a.*Cs1b)./(Cs1a+Cs1b);
Cs1(isnan(Cs1))=0

Cs2a=(eps_pdms.*eps.*(x0+deltaX).*y0)./(z0-deltaZ);
Cs2b=(eps_pdms.*eps.*(x0+deltaX).*y0)./(z0-deltaZ);
Cs2=(Cs2a.*Cs2b)./(Cs2a+Cs2b);

deltaC=(eps_pdms.*eps.*(x0).*y0)./(deltaZ)./2;
figure
plot(deltaZ,Cs1); hold on;
plot(deltaZ,Cs2)
xlabel('Delta-Z (m)')
ylabel('Capacitance (F)')
legend('Cs1','Cs2','Location','best')
title('FSHEAR=+MAX; FNORM=sweep')

% VDIFF mode
Vinp=VDD;
Vinn=VDD;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

figure
subplot(1,2,1)
plot(deltaZ,Vod)
xlabel('Delta-Z (m)')
ylabel('Voutput differential (Voutp-Voutn)')
title('FSHEAR=+MAX; FNORM=sweep')

subplot(1,2,2)
plot(deltaZ,Voutp); hold on
plot(deltaZ,Voutn)
xlabel('Delta-Z (m)')
ylabel('V')
legend('Voutp','Voutn','Location','best')
title('VDIFF MODE (VINP=VINN=VDD)')

Voutp(2)
Voutn(2)
Vod(2)
Vod(end)

% VCM mode
Vinp=VDD;
Vinn=0;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

figure
subplot(1,2,1)
plot(deltaZ,Vod)
xlabel('Delta-Z (m)')
ylabel('Voutput differential (Voutp-Voutn)')
title('FSHEAR=+MAX; FNORM=sweep')

subplot(1,2,2)
plot(deltaZ,Voutp); hold on
plot(deltaZ,Voutn)
xlabel('Delta-Z (m)')
ylabel('V')
legend('Voutp','Voutn','Location','best')
title('VCM MODE (VINP=VDD; VINN=0)')

% VDIFF/VCM mode
Vinp=VDD;
Vinn=VDD;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

Vinp=VDD;
Vinn=0;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vcm=Voutp-Voutn;

Vratio=Vod./Vcm;

figure
subplot(1,2,1)
plot(deltaZ,Vratio)
xlabel('Delta-Z (m)')
ylabel('Voutput ratio (Vod/Vocm)')
title('FSHEAR=+MAX; FNORM=sweep')

subplot(1,2,2)
plot(deltaZ,Vod); hold on
plot(deltaZ,Vcm)
xlabel('Delta-Z (m)')
ylabel('V')
legend('Vod','Vocm','Location','best')
title('VDIFF/VCM RATIOMETRIC MODE')

%% === -MAX CONSTANT SHEAR FORCE; NORMAL FORCE ONLY

x0=50e-6; % overlap at force (shear) = 0
z0=125e-6; % thickness at force (normal) = 0
y0=125e-6; % constant height overlap
deltaZ_max=100e-6; %deltaZ_max=50e-6;
deltaZ_res=deltaZ_max./2^N;
deltaZ=[0:deltaZ_res:deltaZ_max];
deltaX=-50e-6; %max

eps=8.85e-12;
eps_pdms=2.77;
Cs1a=(eps_pdms.*eps.*(x0-deltaX).*y0)./(z0-deltaZ);
Cs1b=(eps_pdms.*eps.*(x0-deltaX).*y0)./(z0-deltaZ);
Cs1=(Cs1a.*Cs1b)./(Cs1a+Cs1b);
Cs1(isnan(Cs1))=0

Cs2a=(eps_pdms.*eps.*(x0+deltaX).*y0)./(z0-deltaZ);
Cs2b=(eps_pdms.*eps.*(x0+deltaX).*y0)./(z0-deltaZ);
Cs2=(Cs2a.*Cs2b)./(Cs2a+Cs2b);
Cs2(isnan(Cs2))=0

deltaC=(eps_pdms.*eps.*(x0).*y0)./(deltaZ)./2;
figure
plot(deltaZ,Cs1); hold on;
plot(deltaZ,Cs2)
xlabel('Delta-Z (m)')
ylabel('Capacitance (F)')
legend('Cs1','Cs2','Location','best')
title('FSHEAR=-MAX; FNORM=sweep')

% VDIFF mode
Vinp=VDD;
Vinn=VDD;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

figure
subplot(1,2,1)
plot(deltaZ,Vod)
xlabel('Delta-Z (m)')
ylabel('Voutput differential (Voutp-Voutn)')
title('FSHEAR=-MAX; FNORM=sweep')

subplot(1,2,2)
plot(deltaZ,Voutp); hold on
plot(deltaZ,Voutn)
xlabel('Delta-Z (m)')
ylabel('V')
legend('Voutp','Voutn','Location','best')
title('VDIFF MODE (VINP=VINN=VDD)')

Voutp(2)
Voutn(2)
Vod(2)
Vod(end)

% VCM mode
Vinp=VDD;
Vinn=0;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

figure
subplot(1,2,1)
plot(deltaZ,Vod)
xlabel('Delta-Z (m)')
ylabel('Voutput differential (Voutp-Voutn)')
title('FSHEAR=-MAX; FNORM=sweep')

subplot(1,2,2)
plot(deltaZ,Voutp); hold on
plot(deltaZ,Voutn)
xlabel('Delta-Z (m)')
ylabel('V')
legend('Voutp','Voutn','Location','best')
title('VCM MODE (VINP=VDD; VINN=0)')

% VDIFF/VCM mode
Vinp=VDD;
Vinn=VDD;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

Vinp=VDD;
Vinn=0;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vcm=Voutp-Voutn;

Vratio=Vod./Vcm;

figure
subplot(1,2,1)
plot(deltaZ,Vratio)
xlabel('Delta-Z (m)')
ylabel('Voutput ratio (Vod/Vocm)')
title('FSHEAR=-MAX; FNORM=sweep')

subplot(1,2,2)
plot(deltaZ,Vod); hold on
plot(deltaZ,Vcm)
xlabel('Delta-Z (m)')
ylabel('V')
legend('Vod','Vocm','Location','best')
title('VDIFF/VCM RATIOMETRIC MODE')
