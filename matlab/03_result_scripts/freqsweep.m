clear variables; close all; clc

load('results/settingsCase1.mat');
load('results/fromCase1.mat');
load('results/fromCase2.mat');
load('results/fromCase3.mat');

% colors of cases 1 2 3

color1 = [179, 179, 179];
color2 = [87, 86, 161];
color3 = [1, 0, 59];

%% Frequency sweep comparison
% VS models - linear
fig1 = figure(1);
fig1.Position = [292 180 759 489];
imp1 = plot(H,FS_ZM1(:,1), 'LineWidth', 1,...
        'Color',[color1(1)/255, color1(2)/255, color1(3)/255]); hold on
imp2 = plot(H,FS_ZM2(:,1), 'LineWidth', 1,...
        'Color',[color2(1)/255, color2(2)/255, color2(3)/255]);
imp3 = plot(H,FS_ZM3(:,1), 'LineWidth', 1,...
        'Color',[color3(1)/255, color3(2)/255, color3(3)/255]);
axis(view)
title('Frequency sweep - VS models - all three cases');
xlabel(xlab);
ylabel(ylab);
legend([imp1,imp2,imp3],'Case 1','Case 2','Case 3');
hold off
clear imp1 imp2 imp3;

% VS models - log
fig2 = figure(2);
fig2.Position = [292 180 759 489];
imp1 = semilogy(H,FS_ZM1(:,1), 'LineWidth', 1,...
        'Color',[color1(1)/255, color1(2)/255, color1(3)/255]); hold on
imp2 = semilogy(H,FS_ZM2(:,1), 'LineWidth', 1,...
        'Color',[color2(1)/255, color2(2)/255, color2(3)/255]);
imp3 = semilogy(H,FS_ZM3(:,1), 'LineWidth', 1,...
        'Color',[color3(1)/255, color3(2)/255, color3(3)/255]);
title('Frequency sweep - VS models - all three cases');
xlabel(xlab);
ylabel(strcat(ylab,' (log axis)'));
legend([imp1,imp2,imp3],'Case 1','Case 2','Case 3');
hold off
clear imp1 imp2 imp3;

% Z(s) models - linear
fig3 = figure(3);
fig3.Position = [292 180 759 489];
imp1 = plot(H,FS_ZM1(:,2), 'LineWidth', 1,...
        'Color',[color1(1)/255, color1(2)/255, color1(3)/255]); hold on
imp2 = plot(H,FS_ZM2(:,2), 'LineWidth', 1,...
        'Color',[color2(1)/255, color2(2)/255, color2(3)/255]);
imp3 = plot(H,FS_ZM3(:,2), 'LineWidth', 1,...
        'Color',[color3(1)/255, color3(2)/255, color3(3)/255]);
axis(view)
title('Frequency sweep - Z(s) models - all three cases');
xlabel(xlab);
ylabel(ylab);
legend([imp1,imp2,imp3],'Case 1','Case 2','Case 3');
hold off
clear imp1 imp2 imp3;

% Z(s) models - log
fig4 = figure(4);
fig4.Position = [292 180 759 489];
imp1 = semilogy(H,FS_ZM1(:,2), 'LineWidth', 1,...
        'Color',[color1(1)/255, color1(2)/255, color1(3)/255]); hold on
imp2 = semilogy(H,FS_ZM2(:,2), 'LineWidth', 1,...
        'Color',[color2(1)/255, color2(2)/255, color2(3)/255]);
imp3 = semilogy(H,FS_ZM3(:,2), 'LineWidth', 1,...
        'Color',[color3(1)/255, color3(2)/255, color3(3)/255]);;
title('Frequency sweep - Z(s) models - all three cases');
xlabel(xlab);
ylabel(strcat(ylab,' (log axis)'));
legend([imp1,imp2,imp3],'Case 1','Case 2','Case 3');
hold off
clear imp1 imp2 imp3;

%% HMA comparison
% VS models
fig5 = figure(5);
fig5.Position = [244 115 847 576];
hma1 = plot(H,HMA_ZmaxM1(:,1), 'LineWidth', 1,...
        'Color',[color1(1)/255, color1(2)/255, color1(3)/255]); hold on
hma2 = plot(H,HMA_ZmaxM2(:,1), 'LineWidth', 1,...
        'Color',[color2(1)/255, color2(2)/255, color2(3)/255]);
hma3 = plot(H,HMA_ZmaxM3(:,1), 'LineWidth', 1,...
        'Color',[color3(1)/255, color3(2)/255, color3(3)/255]);
axis(view)
title('HMA - VS models - all three cases');    
xlabel(xlab);
ylabel(ylab);
legend([hma1,hma2,hma3],'Case 1','Case 2','Case 3');
hold off
clear hma1 hma2 hma3;

% VS models
fig6 = figure(6);
fig6.Position = [244 115 847 576];
hma1 = plot(H,HMA_ZmaxM1(:,2), 'LineWidth', 1,...
        'Color',[color1(1)/255, color1(2)/255, color1(3)/255]); hold on
hma2 = plot(H,HMA_ZmaxM2(:,2), 'LineWidth', 1,...
        'Color',[color2(1)/255, color2(2)/255, color2(3)/255]);
hma3 = plot(H,HMA_ZmaxM3(:,2), 'LineWidth', 1,...
        'Color',[color3(1)/255, color3(2)/255, color3(3)/255]);
axis(view)
title('HMA - Z(s) models - all three cases');    
xlabel(xlab);
ylabel(ylab);
legend([hma1,hma2,hma3],'Case 1','Case 2','Case 3');
hold off
clear hma1 hma2 hma3;

