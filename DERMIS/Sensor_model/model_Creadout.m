clear all
close all

C0=0.6125e-15; % no normal force
N=5;
deltaC=[0:C0/2^N:C0];
Cf=10e-15;
Cs1=C0-deltaC; Cs2=C0+deltaC;
VDD=3.3;
Vocm=VDD/2;

figure
plot(deltaC,Cs1); hold on;
plot(deltaC,Cs2)

% VDIFF mode
Vinp=VDD;
Vinn=VDD;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

figure
subplot(1,2,1)
plot(deltaC,Vod)
subplot(1,2,2)
plot(deltaC,Voutp); hold on
plot(deltaC,Voutn)

Voutp(2)
Voutn(2)
Vod(2)

% VCM mode
Vinp=VDD;
Vinn=0;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

figure
subplot(1,2,1)
plot(deltaC,Vod)
subplot(1,2,2)
plot(deltaC,Voutp); hold on
plot(deltaC,Voutn)

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
plot(deltaC,Vratio)
subplot(1,2,2)
plot(deltaC,Vod); hold on
plot(deltaC,Vcm)