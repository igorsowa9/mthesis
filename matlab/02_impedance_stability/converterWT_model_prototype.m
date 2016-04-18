clear variables; close all; clc

%% WT converter
Sr = 5e3;
V0 = 400; % dc bus 
V1 = 170; % ph. volt. ampl. 
I1 = Sr/(3*V1);

L = 3.4e-3; % ph. ind.

% Hi(s) current control compensator
Kp_i = 0.047;
Ki_i = 148;

% Hp(s) PLL compensator
Kp_p = 0.269;
Ki_p = 16.42;

f1 = 50;
w1 = 2*pi*f1;
res = 0.01;
H = 0.2:res:200; % 10Hz - 10kHz
Hf = H*f1;

Hi = @(s) Kp_i+Ki_i/s;
Hp = @(s) (Kp_p+Ki_i/s)/s;

ZpM = zeros(1,length(H));
ZnM = zeros(1,length(H));
for hh=1:length(H)
    h=H(hh);
    s = h*1i*w1;
    
    Tpll_subt1 = V1*Hp(s-1i*w1)/(2*(1+V1*Hp(s-1i*w1)));
    ZpM(hh) = (Hi(s-1i*w1)*V0 + (s-1i*w1)*L)/...
        (1-Tpll_subt1*(1+Hi(s-1i*w1)*I1*V0/V1));
    
    Tpll_add1 = V1*Hp(s+1i*w1)/(2*(1+V1*Hp(s+1i*w1)));
    ZnM(hh) = (Hi(s+1i*w1)*V0 + (s+1i*w1)*L)/...
        (1-Tpll_add1*(1+Hi(s+1i*w1)*I1*V0/V1));
end
figure(1)
semilogx(Hf,mag2db(abs(ZpM))); hold on
semilogx(Hf,mag2db(abs(ZnM)), 'r'); 
hold off

figure(2)
semilogx(Hf,rad2deg(angle(ZpM))); hold on
semilogx(Hf,rad2deg(angle(ZnM)), 'r'); 
hold off



%%%%