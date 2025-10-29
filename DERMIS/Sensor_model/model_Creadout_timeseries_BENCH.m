clear all
close all

%% ==== NO NORMAL FORCE =====
N = 5;
Cf=40e-15;
Coffset1=0%0.85e-12 %1e-15;
Coffset2=0%0.85e-12 %1e-15;

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
eps_pdms = 2.77; %factor 100x larger capacitance

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

Ccancel1=0%Cs2%0%Cs1;
Ccancel2=0%Cs1%0%Cs2;
beta1=Cf./(Cf+Cs2+Ccancel1); beta2=Cf./(Cf+Cs1+Ccancel2);
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
%plot(t,1./(beta1+beta2) .* (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*(1/Ao+beta1)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
legend('VTX component','VOCM component','Location','best')

subplot(1,2,2)
plot(t,1./(beta1+beta2) .* (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2))  ) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
hold on;
plot(t,1./(beta1+beta2) .* ( 2*Vocm.*(1/Ao+beta2)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
plot(t,Voutn)
%plot(t,1./(beta1+beta2) .* (1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2))  + 2*Vocm.*(1/Ao+beta2)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
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



%% ==== NO NORMAL FORCE (BUT WITH BENCH CAPACITIVE INPUT) =====

%close all
clear all


N = 5;
Cf=40e-15;
Coffset1=1.7e-12%0.85e-12 %1e-15;
Coffset2=1.7e-12%0.85e-12 %1e-15;

% Cs1=C0-deltaC; Cs2=C0+deltaC;
VDD=3 %3.3;
Vocm=1.5 %VDD/2; % this is the source of nonlinearity for VCM measurements
Ao=55 %100;

fin = 100; % Define the frequency of the sinusoid (1/fin)
deltaX_max = 50e-6;
deltaX_res = deltaX_max / (2^N);

% Create a sinusoidal waveform for deltaX
t = linspace(0, 1, 2^12); % Time vector for one period
C_amp = (4.33e-12 - Coffset1);
Cs1 = Coffset1 + (C_amp/2) + (C_amp/2) * sin(2 * pi * fin * t); % Sinusoidal waveform
Cs2 = Coffset1 + (C_amp/2) - (C_amp/2) * sin(2 * pi * fin * t); % Sinusoidal waveform
% deltaX = deltaX_max * sin(2 * pi * fin * t); % Sinusoidal waveform
% 
% % Repeat the waveform to cover the range from -deltaX_max to deltaX_max
% % deltaX = interp1(linspace(0, 1, length(deltaX)), deltaX, linspace(0, 1, 2^N + 1));
% % deltaX = repmat(deltaX, 1, floor((2 * deltaX_max / fin)));
% 
% x0 = 50e-6;
% z0 = 125e-6;
% y0 = 125e-6;
% deltaZ = 100e-6 %0;
% 
% eps = 8.85e-12;
% eps_pdms = 2.77*500; %factor 300x larger capacitance
% 
% % Calculation for Cs1 and Cs2 using deltaX as a sinusoidal waveform
% Cs1a = (eps_pdms * eps * (x0 + deltaX) * y0) / (z0 - deltaZ);
% Cs1b = (eps_pdms * eps * (x0 + deltaX) * y0) / (z0 - deltaZ);
% %Cs1a = 100.*Cs1a; Cs1b = 100.*Cs1b;
% Cs1 = Coffset1 + (Cs1a .* Cs1b) ./ (Cs1a + Cs1b);
% %Cs1(isnan(Cs1)) = 0;
% 
% Cs2a = (eps_pdms * eps * (x0 - deltaX) * y0) / (z0 - deltaZ);
% Cs2b = (eps_pdms * eps * (x0 - deltaX) * y0) / (z0 - deltaZ);
% %Cs2a = 100.*Cs2a; Cs2b = 100.*Cs2b;
% Cs2 = Coffset2 +(Cs2a .* Cs2b) ./ (Cs2a + Cs2b);
% Cs2(isnan(Cs2)) = 0;

% figure
% plot(t,deltaX)
% legend('deltaX','Location','best')

figure
% plot(t,Cs1a)
% hold on
plot(t,Cs1)
% legend('Cs1a (single)','Cs1 (parallel)','Location','best')

figure
% plot(t,Cs2a)
% hold on
plot(t,Cs2)
% legend('Cs2a (single)','Cs2 (parallel)','Location','best')

    
% VDIFF mode
Vinp=VDD; %0 IF CT; VDD IF DT
Vinn=VDD; %0 IF CT; VDD IF DT

Ccancel1=0%Cs2%0%Cs1;
Ccancel2=0%Cs1%0%Cs2;
beta1=Cf./(Cf+Cs2+Ccancel2); beta2=Cf./(Cf+Cs1+Ccancel1);
beta1_in= Cs2./(Cf+Cs2+Ccancel2); beta2_in= Cs1./(Cf+Cs1+Ccancel1);
%Voutp = 1./(beta1+beta2) .* ((Vinp).*(beta1_in) - (Vinn).*(beta2_in) + 2*Vocm.*(1/Ao+beta1)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2));
Voutp = 1./(beta1+beta2) .* ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*(1/Ao+beta1)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2));
%Voutn= 1./(beta1+beta2) .* (-1*((Vinp).*(beta1_in) - (Vinn).*(beta2_in))  + 2*Vocm.*(1/Ao+beta2)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2));
Voutn= 1./(beta1+beta2) .* (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2))  + 2*Vocm.*(1/Ao+beta2)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2));
Vod=Voutp-Voutn;
%Vod can be approximated as 2*(Vinp-Vinn)*(Cs/Cf) if Cs=Cs1=Cs2
%Voutp can be approx. as (Vinp-Vinn)*(Cs+delC)/Cf
%Voutn can be approx. as -(Vinp-Vinn)*(Cs-delC)/Cf

figure
subplot(1,2,1)
plot(t,1./(beta1+beta2) .* ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) ) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
%plot(t,1./(beta1+beta2) .* ((Vinp).*(beta1_in) - (Vinn).*(beta2_in) ) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
hold on;
% DEPRACATED: 1/Ao more comparable to beta if beta is small (Cs1 >> Cf)
% Voutp = (2*Vocm/2.xx) * beta1 / (beta1+beta2) ->> reason why output is working even with pF range
% input
% plot(t,Vocm.*(beta1./(beta1+beta2)))
plot(t,1./(beta1+beta2) .* (2*Vocm.*(1/Ao+beta1)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
% Voutp is  [Vocm - A*sin(2pift)] + 2A*sin(2pift) = Vocm - A*sin(2pift) ->
% inverse of Vocm component
plot(t,Voutp)
%plot(t,1./(beta1+beta2) .* (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*(1/Ao+beta1)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
legend('VTX component','VOCM component','Location','best')

subplot(1,2,2)
%plot(t,1./(beta1+beta2) .* (-1*((Vinp).*(beta1_in) - (Vinn).*(beta2_in))  ) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
plot(t,1./(beta1+beta2) .* (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2))  ) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
hold on;
plot(t,1./(beta1+beta2) .* ( 2*Vocm.*(1/Ao+beta2)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
plot(t,Voutn)
%plot(t,1./(beta1+beta2) .* (1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2))  + 2*Vocm.*(1/Ao+beta2)) ./ (1 + 2./(Ao.*beta1+Ao.*beta2)))
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

figure
plot(t,Vocm.*(beta1./(beta1+beta2)))
figure
plot(t,Cs1./(Cs1+Cs2))


%% what happens to VP and VN
Vp = Vinp .* (Cs1 ./(Cs1+Ccancel1+Cf)) + Voutn .* (Cf ./(Cs1+Ccancel1+Cf))
figure 
plot(t,Vp)