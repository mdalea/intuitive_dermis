%texsignal = csvread('7p20mvpp-supplyext.csv');

texsignal = csvread('/users/micas/malea/Documents/ArduinoRead_2/ArduinoRead/7p20mvpp-supplyext_nospks-ch0-bw10khz.csv');

%afe_filt = tf([1e3 0],[1 1000.5 500]);
%[texsignal_filt,~] = lsim(afe_filt,texsignal(:,2),texsignal(:,1));
period = 0.1;
origperiod = 0.4;
timenorm = origperiod/period;

figure
%plot(timenorm.*(texsignal(:,1)-texsignal(1,1))-0.0625,texsignal(:,2))
plot(texsignal(:,1),texsignal(:,2))
%xlim([0 0.4])
xlabel('Time (s)')
ylabel('Texture signal')
micasplot_large




