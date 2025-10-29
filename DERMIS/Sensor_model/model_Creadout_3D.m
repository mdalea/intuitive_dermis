% remember this only works for AC VTX, results are AC magnitudes
%Cs1 and Cs2 are reversed but who cares

% 3D plot vs deltaX and deltaZ of the sensor

clear all
close all

%% ==== SENSOR/READOUT PARAMS =====
% C0=0.6125e-15; % no normal force
N=5;
% deltaC=[0:C0/2^N:C0];
Cf=10e-15;
% Cs1=C0-deltaC; Cs2=C0+deltaC;
VDD=3.3;
Vocm=VDD/2; % this is the source of nonlinearity for VCM measurements


%% === NO SHEAR FORCE; NORMAL FORCE ONLY

x0=50e-6; % overlap at force (shear) = 0
z0=125e-6; % thickness at force (normal) = 0
y0=125e-6; % constant height overlap
deltaZ_max=100e-6; %deltaZ_max=50e-6;
deltaZ_res=deltaZ_max./2^N./2; % factor 2 to have same length as deltaX
deltaZ=[0:deltaZ_res:deltaZ_max];
deltaX_max=50e-6;
deltaX_res=deltaX_max./2^N;
deltaX=[-deltaX_max:deltaX_res:deltaX_max];

% Create grid for surf plot
[DeltaX, DeltaZ] = meshgrid(deltaX, deltaZ);

eps=8.85e-12;
eps_pdms=2.77;%*1000;
Cs1a=(eps_pdms.*eps.*(x0-DeltaX).*y0)./(z0-DeltaZ);
Cs1b=(eps_pdms.*eps.*(x0-DeltaX).*y0)./(z0-DeltaZ);
Cs1=(Cs1a.*Cs1b)./(Cs1a+Cs1b);
Cs1(isnan(Cs1))=0

Cs2a=(eps_pdms.*eps.*(x0+DeltaX).*y0)./(z0-DeltaZ);
Cs2b=(eps_pdms.*eps.*(x0+DeltaX).*y0)./(z0-DeltaZ);
Cs2=(Cs2a.*Cs2b)./(Cs2a+Cs2b);

deltaC=(eps_pdms.*eps.*(x0).*y0)./(DeltaZ)./2;
% Reshape Cs1 to match the size of the grid
% Cs1_matrix = reshape(Cs1, numel(deltaZ), numel(deltaX));

% Plotting surface
figure;
surf(DeltaX, DeltaZ, Cs1);
hold on;
surf(DeltaX, DeltaZ, Cs2);

%plot3(deltaX, deltaZ, Cs1, 'b-', 'LineWidth', 1.5);

% plot(deltaZ,Cs1); hold on;
% plot(deltaZ,Cs2)
xlabel('Delta-X (m)')
ylabel('Delta-Z (m)')
zlabel('Capacitance (F)')
legend('Cs1','Cs2','Location','best')
%title('FSHEAR=0; FNORM=sweep')

% VDIFF mode
Vinp=VDD;
Vinn=VDD;

beta1=Cf./(Cf+Cs2); beta2=Cf./(Cf+Cs1);
Voutp= ((Vinp).*(1-beta1) - (Vinn).*(1-beta2) + 2*Vocm.*beta1) ./ (beta1+beta2);
Voutn= (-1*((Vinp).*(1-beta1) - (Vinn).*(1-beta2)) + 2*Vocm.*beta2) ./ (beta1+beta2);
Vod=Voutp-Voutn;

figure
subplot(1,2,1)
%plot(deltaZ,Vod)
surf(DeltaX, DeltaZ, Vod);
xlabel('Delta-X (m)')
ylabel('Delta-Z (m)')
zlabel('Voutput differential (Voutp-Voutn)')
%title('FSHEAR=0; FNORM=sweep')

subplot(1,2,2)
% plot(deltaZ,Voutp); hold on
% plot(deltaZ,Voutn)
surf(DeltaX, DeltaZ, Voutp); hold on
surf(DeltaX, DeltaZ, Voutn);
xlabel('Delta-X (m)')
ylabel('Delta-Z (m)')
zlabel('V')
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
surf(DeltaX, DeltaZ, Vod); 
xlabel('Delta-X (m)')
ylabel('Delta-Z (m)')
zlabel('Voutput differential (Voutp-Voutn)')
%title('FSHEAR=0; FNORM=sweep')

subplot(1,2,2)
surf(DeltaX, DeltaZ, Voutp); hold on
surf(DeltaX, DeltaZ, Voutn);
xlabel('Delta-X (m)')
ylabel('Delta-Z (m)')
zlabel('V')
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
%plot(deltaZ,Vratio)
surf(DeltaX, DeltaZ, Vratio); 
xlabel('Delta-X (m)')
ylabel('Delta-Z (m)')
zlabel('Voutput ratio (Vod/Vocm)')
%title('FSHEAR=0; FNORM=sweep')

subplot(1,2,2)
% plot(deltaZ,Vod); hold on
% plot(deltaZ,Vcm)
surf(DeltaX, DeltaZ, Vod); hold on
surf(DeltaX, DeltaZ, Vcm);
xlabel('Delta-X (m)')
ylabel('Delta-Z (m)')
zlabel('V')
legend('Vod','Vocm','Location','best')
title('VDIFF/VCM RATIOMETRIC MODE')

%%

figure
plot(Cs1(1,:),Cs1(1,:)./Vod)