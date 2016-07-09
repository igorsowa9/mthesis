clear all
close all
clc

L1 = 3.8e-3;
L2 = 1.9e-3;
C = 79.44e-6;
R = 1.3353;

s = tf('s');
Hlcl = (s*C*R+1)/(s^2*C*L1+s*C*R+1);
bode(Hlcl);
grid on