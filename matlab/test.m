clear all
close all
clc

plot(-2:5, (-2:5).^2-1)

%# vertical line
hx = graph2d.constantline(0, 'LineStyle',':', 'Color',[.7 .7 .7]);
changedependvar(hx,'x');

%# horizontal line
hy = graph2d.constantline(0, 'Color',[.7 .7 .7]);
changedependvar(hy,'y');

break



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