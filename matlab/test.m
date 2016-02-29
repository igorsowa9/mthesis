clear all
close all
clc

f = 50;
tmax = 0.1;
t = 0:0.0001:tmax;
A = 1;
w = 2*pi*f;

V1 = A*1*sin(w*t);
V3 = A*0.11*sin(w*3*t);
V5 = A*0.08*sin(w*5*t);
V9 = A*0.07*sin(w*9*t);

V=V1+V3+V5+V9;

figure(1)
plot(t,V);

figure(2)
plot(10*(0:length(t)-1),abs(fft(V))/f*tmax);
%axis([0 50 0 15000]);