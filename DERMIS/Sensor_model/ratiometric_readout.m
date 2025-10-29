%% understanding ratiometric shit
delta_d=[1:1:9999]; % change in displacement
d0=10000; % constant displacement

x=delta_d./d0 % needs to be <1
CS1=C0./(1-x)
CS2=C0./(1+x)

figure;plot(x,CS1); hold on;
plot(x,CS2)
title('CS1 CS2')

%non-ratiometric differential
figure;plot(x,CS1-CS2)
title('CS1-CS2')

%ratiometric differential
figure;plot(x,(CS1-CS2)./(CS1+CS2))
title('ratiometric')