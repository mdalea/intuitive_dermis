clear all
close all

C0=6e-15;
Cdac=[1.7e-12:0.38e-12:193.84e-12]-1.7e-12;
Cseries=(C0.*Cdac)./(C0+Cdac)

figure
plot(Cdac,Cseries,'-o')