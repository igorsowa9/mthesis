clear vars
clc
close all


H = 0.01:0.01:30;
f = H*50;
w = 2*pi*f;

R = 4;
L = 12.5e-3;
C = 10e-6;

%% series resonance
Zs = R+1i*(w*L-1./(w*C));

figure(1)
subplot(2,1,1)
semilogy(H,abs(Zs), 'LineWidth', 2);
title('Series Resonance');
ylabel('Magnitude [Ohm]');
grid on

subplot(2,1,2)
plot(H,rad2deg(angle(Zs)), 'LineWidth', 2);
xlabel('Order');
ylabel('Angle [deg]');
grid on

f_series = (1/sqrt(L*C))/(2*pi)

%% parallel resonance

Yp = 1/R + 1./(1i*w*L) + 1i*w*C;
Zp = 1./Yp;

figure(2)
subplot(2,1,1)
semilogy(H,abs(Zp), 'LineWidth', 2);
title('Parallel Resonance');
ylabel('Magnitude [Ohm]');
grid on

subplot(2,1,2)
plot(H,rad2deg(angle(Zp)), 'LineWidth', 2);
xlabel('Order');
ylabel('Angle [deg]');
grid on

f_parallel = (1/sqrt(L*C))/(2*pi)

%% tank circuit

% Yt = R./(R^2+w.^2*L^2) - 1i*(w*L./(R^2+w.^2*L^2) - w*C);
Yt = 1i*w*C + 1./(R+1i*w*L);
Zt = 1./Yt;

figure(3)
subplot(2,1,1)
semilogy(H,abs(Zt), 'LineWidth', 2);
title('Parallel Tank Circuit Resonance');
ylabel('Magnitude [Ohm]');
grid on

subplot(2,1,2)
plot(H,rad2deg(angle(Zt)), 'LineWidth', 2);
xlabel('Order');
ylabel('Angle [deg]');
grid on

f_tank = sqrt(1/(L*C) - R^2/L^2)/(2*pi)


