clear all
%close all

%% ==== NO NORMAL FORCE =====
N = 5;
Cf=40e-15;
Coffset1=0 %1e-15;
Coffset2=0 %1e-15;
% Cs1=C0-deltaC; Cs2=C0+deltaC;
VDD=3 %3.3;
Vocm=1.5 %VDD/2; % this is the source of nonlinearity for VCM measurements
Ao=40 %100;

fin = 100; % Define the frequency of the sinusoid (1/fin)
deltaX_max = 50e-6;
deltaX_res = deltaX_max / (2^N);

% Create a sinusoidal waveform for deltaX
t = linspace(0, 1, 2^12); % Time vector for one period
deltaX = deltaX_max * sin(2 * pi * fin * t); % Sinusoidal waveform

% Repeat the waveform to cover the range from -deltaX_max to deltaX_max
% deltaX = interp1(linspace(0, 1, length(deltaX)), deltaX, linspace(0, 1, 2^N + 1));
% deltaX = repmat(deltaX, 1, floor((2 * deltaX_max / fin)));

x0 = 50e-6;
z0 = 125e-6;
y0 = 125e-6;
deltaZ = 100e-6 %0;

eps = 8.85e-12;
eps_pdms = 2.77;

% Calculation for Cs1 and Cs2 using deltaX as a sinusoidal waveform
Cs1a = (eps_pdms * eps * (x0 + deltaX) * y0) / (z0 - deltaZ);
Cs1b = (eps_pdms * eps * (x0 + deltaX) * y0) / (z0 - deltaZ);
%Cs1a = 100.*Cs1a; Cs1b = 100.*Cs1b;
Cs1 = Coffset1 + (Cs1a .* Cs1b) ./ (Cs1a + Cs1b);
Cs1(isnan(Cs1)) = 0;

Cs2a = (eps_pdms * eps * (x0 - deltaX) * y0) / (z0 - deltaZ);
Cs2b = (eps_pdms * eps * (x0 - deltaX) * y0) / (z0 - deltaZ);
%Cs2a = 100.*Cs2a; Cs2b = 100.*Cs2b;
Cs2 = Coffset2 +(Cs2a .* Cs2b) ./ (Cs2a + Cs2b);
Cs2(isnan(Cs2)) = 0;

figure
plot(t,deltaX)
legend('deltaX','Location','best')

figure
plot(t,Cs1a)
hold on
plot(t,Cs1)
legend('Cs1a (single)','Cs1 (parallel)','Location','best')

figure
plot(t,Cs2a)
hold on
plot(t,Cs2)
legend('Cs2a (single)','Cs2 (parallel)','Location','best')


% VDIFF mode
Vinp=VDD; %0 IF CT; VDD IF DT
Vinn=VDD; %0 IF CT; VDD IF DT

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
%Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutp = 1./(beta1+beta2) .* ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*(1/Ao+beta1)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2));
%Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Voutn= 1./(beta1+beta2) .* (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2))  + 2*Vocm.*(1/Ao+beta2)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2));
Vod=Voutp-Voutn;
%Vod can be approximated as 2*(Vinp-Vinn)*(Cs/Cf) if Cs=Cs1=Cs2
%Voutp can be approx. as (Vinp-Vinn)*(Cs+delC)/Cf
%Voutn can be approx. as -(Vinp-Vinn)*(Cs-delC)/Cf

figure
subplot(1,2,1)
plot(t,1./(beta1+beta2) .* ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) ) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
hold on;
plot(t,1./(beta1+beta2) .* (2*Vocm.*(1/Ao+beta1)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
plot(t,Voutp)
legend('VTX component','VOCM component','Location','best')

subplot(1,2,2)
plot(t,1./(beta1+beta2) .* (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2))  ) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
hold on;
plot(t,1./(beta1+beta2) .* ( 2*Vocm.*(1/Ao+beta2)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
plot(t,Voutn)
legend('VTX component','VOCM component','Location','best')

figure

subplot(1,2,1)
plot(t,Vod)
xlabel('time (s)')
ylabel('Voutput differential (Voutp-Voutn)')
title('FSHEAR=sweep; FNORM=0')

subplot(1,2,2)
plot(t,Voutp); hold on
plot(t,Voutn)
xlabel('time (s)')
ylabel('V')
legend('Voutp','Voutn','Location','best')
title('VDIFF MODE (VINP=VINN=VDD)')

% VCM mode
Vinp=VDD %0 %VDD;
Vinn=0 %VDD %0;


beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
%Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutp = 1./(beta1+beta2) .* ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*(1/Ao+beta1)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2));
%Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Voutn= 1./(beta1+beta2) .* (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2))  + 2*Vocm.*(1/Ao+beta2)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2));
Vod=Voutp-Voutn;

figure

subplot(1,2,1)
plot(t,Vod)
xlabel('time (s)')
ylabel('Voutput differential (Voutp-Voutn)')
title('FSHEAR=sweep; FNORM=0')

subplot(1,2,2)
plot(t,Voutp); hold on
plot(t,Voutn)
xlabel('time (s)')
ylabel('V')
legend('Voutp','Voutn','Location','best')
title('VCM MODE (VINP=VINN=out of phase)')

