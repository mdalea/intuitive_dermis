% Derive numerically the Cs1 - Cs2 (deltaC) transfer function to VOD,
% including capacitance-referred noise RMS /noise floor

% remember this only works for AC VTX, results are AC magnitudes
%Cs1 and Cs2 are reversed but who cares

clear all
close all


%% ==== SHEAR FORCE SWEEP, NO NORMAL FORCE =====
clear all
close all

% C0=0.6125e-15; % no normal force
N=5;
% deltaC=[0:C0/2^N:C0];
Cf=10e-15;
% Cs1=C0-deltaC; Cs2=C0+deltaC;
VDD=3;
Vocm=VDD/2; % this is the source of nonlinearity for VCM measurements
% fhi=1e3; flo=1;
fsig_bw=1e3;
fhi=(pi/2)*fsig_bw*5*5; % fbw_ota = 5 * ftx; ftx = 5 * fsig_bw
flo=0;

x0=50e-6; % overlap at force (shear) = 0
z0=125e-6; % thickness at force (normal) = 0
y0=125e-6; % constant height overlap
deltaX_max=50e-6;
deltaX_res=deltaX_max./2^N;
%deltaX=[-deltaX_max:deltaX_res:deltaX_max];
deltaX=[deltaX_res:deltaX_res:deltaX_max];
deltaZ=0;
Vir_noise_rms = linspace(10e-9, 85e-6, numel(deltaX))   %100uV worst case -> to make AFE noise << mismatch in comparator

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

figure
%plot(Vod,sqrt((1e-3)^2.*  (1./beta1).^2 .* ((Cs1-Cs2)./Vod).^2))
plot(Vod,(Cs1-Cs2)./Vod)
title('Transfer function (Cs1-Cs2)/Vod - VDIFF mode')


% DOUBLE-CHECK EQUATIONS HERE. MADE IN A HURRY
figure
Vout_noise_rms = sqrt(Vir_noise_rms.^2 .* max((1./beta1)).^2);
Cs1_minus_Cs2_referred_noise_rms = sqrt(Vout_noise_rms.^2 .* ((Cs1-Cs2)./Vod).^2);
% Vout_noise_floor = Vout_noise_rms  .* sqrt(fhi - flo);
Vout_noise_floor = Vout_noise_rms  ./ sqrt(fhi - flo);
plot(Vout_noise_floor, Cs1_minus_Cs2_referred_noise_rms)
ylabel('Cs1 - Cs2 referred noise RMS')
xlabel('Vinput referred noise floor (V/sqrt(Hz)) ')
title('VDIFF mode')

figure
plot(Vout_noise_floor.^2, Cs1_minus_Cs2_referred_noise_rms)
ylabel('Cs1 - Cs2 referred noise RMS')
xlabel('Vinput referred noise floor (V^2/Hz) ')
title('VDIFF mode')

figure
% gm = (2*4 * (2/3) * 1.38e-23 * 300) ./ (Vout_noise_floor.^2); % differential
gm = (2*4 * (2/3) * 1.38e-23 * 300) ./ ((Vir_noise_rms.^2)/(fhi-flo)); % differential
gm(end)
plot(gm, Cs1_minus_Cs2_referred_noise_rms)
ylabel('Cs1 - Cs2 referred noise RMS')
xlabel('Required Gm (assuming input pair noise only) ')
title('VDIFF mode')

figure
Vout_noise_rms = sqrt(Vir_noise_rms.^2 .* (max(1./beta1)).^2);
Cs1_minus_Cs2_referred_noise_rms = sqrt(Vout_noise_rms.^2 .* ((Cs1-Cs2)./Vod).^2);
plot(Vir_noise_rms, Cs1_minus_Cs2_referred_noise_rms)
xlabel('Vinput referred noise RMS')
ylabel('Cs1 - Cs2 referred noise RMS')
title('VDIFF mode')

figure
% fhi=1e3; flo=1;
% Cs1_minus_Cs2_referred_noise_floor = Cs1_minus_Cs2_referred_noise_rms .* sqrt(fhi - flo);
Cs1_minus_Cs2_referred_noise_floor = Cs1_minus_Cs2_referred_noise_rms ./ sqrt(fhi - flo);
plot(Vir_noise_rms, Cs1_minus_Cs2_referred_noise_floor)
xlabel('Vinput referred noise RMS')
ylabel('Cs1 - Cs2 referred noise floor (F/sqrt(Hz)')
title('VDIFF mode')

figure
RC_kTC_noise_rms = linspace(1e-6, 100e-6, numel(deltaX));
Cs1_minus_Cs2_referred_noise_rms_kTC = sqrt(RC_kTC_noise_rms.^2 .* ((Cs1-Cs2)./Vod).^2);
plot(RC_kTC_noise_rms, Cs1_minus_Cs2_referred_noise_rms_kTC)
xlabel('RC filter kTC noise RMS')
ylabel('Cs1 - Cs2 referred noise RMS')
title('VDIFF mode')

%% VCM mode
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


figure
%plot(Vod,sqrt((1e-3)^2.*  (1./beta1).^2 .* ((Cs1-Cs2)./Vod).^2))
plot(Vod,(Cs1-Cs2)./Vod)
title('Transfer function (Cs1-Cs2)/Vod - VCM mode')

figure
Vout_noise_rms = sqrt(Vir_noise_rms.^2 .* (1./beta1).^2);
Cs1_minus_Cs2_referred_noise_rms = sqrt(Vout_noise_rms.^2 .* ((Cs1-Cs2)./Vod).^2)
plot(Vir_noise_rms, Cs1_minus_Cs2_referred_noise_rms)
xlabel('Vinput referred noise RMS')
ylabel('Cs1 - Cs2 referred noise RMS')
title('VCM mode')

figure
% fhi=1e3; flo=1;
% Cs1_minus_Cs2_referred_noise_floor = Cs1_minus_Cs2_referred_noise_rms .* sqrt(fhi - flo)
Cs1_minus_Cs2_referred_noise_floor = Cs1_minus_Cs2_referred_noise_rms ./ sqrt(fhi - flo)
plot(Vir_noise_rms, Cs1_minus_Cs2_referred_noise_floor)
xlabel('Vinput referred noise RMS')
ylabel('Cs1 - Cs2 referred noise floor (F/sqrt(Hz)')
title('VCM mode')


%% VDIFF/VCM mode
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

figure
%plot(Vod,sqrt((1e-3)^2.*  (1./beta1).^2 .* ((Cs1-Cs2)./Vod).^2))
plot(Vod./Vcm,(Cs1-Cs2)./(Vod./Vcm))
title('Transfer function (Cs1-Cs2)/Vod - VCM mode')

figure
Vout_noise_rms = sqrt(Vir_noise_rms.^2 .* (1./beta1).^2);
Cs1_minus_Cs2_referred_noise_rms = sqrt(Vout_noise_rms.^2 .* ((Cs1-Cs2)./(Vod./Vcm)).^2)
plot(Vir_noise_rms, Cs1_minus_Cs2_referred_noise_rms)
xlabel('Vinput referred noise RMS')
ylabel('Cs1 - Cs2 referred noise RMS')
title('VDIFF/VCM mode')

%% COMPARATOR
fullscale_pp=2*230e-3;
LSB=fullscale_pp./2^7;
maxSR=2*pi*(fullscale_pp/2)*fsig_bw;
tloop=LSB./maxSR;
fbw_comp=3/(2*pi*tloop);
fhi_comp=(pi/2)*fbw_comp; % fbw_ota = 5 * ftx; ftx = 5 * fsig_bw
flo_comp=0;

Vir_noise_comp_rms = linspace(10e-9, 0.1*LSB, numel(deltaX))   %100uV worst case -> to make AFE noise << mismatch in comparator

gm_comp = (4 * (2/3) * 1.38e-23 * 300) ./ ((Vir_noise_comp_rms.^2)/(fhi_comp-flo_comp)); % differential
gm_comp(end)