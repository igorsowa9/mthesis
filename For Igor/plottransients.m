%% plot waveform transients (AC currents, DC voltage, DC current)

[fname, pname] = uigetfile('*.wvf', 'Select WVF files');
filename = [pname fname];

[iL1,t] = wvfread(filename, 1, 2);
iL2 = wvfread(filename, 1, 4);
iL3 = wvfread(filename, 2, 2);
uDC = wvfread(filename, 2, 3);
iDC = wvfread(filename, 2, 4);

info = hdrread(filename);
timestring = info.Group1.Trace1.Time{1,1};

% scaling correction
uDC = 2*uDC;

% plots 
figure('Name', ['Transients at ', timestring])
clf
plotbrowser('on')

%%% plot with 3 subplots

aX1 = subplot(3,1,1);
plot(t, [iL1 iL2 iL3])
xlabel('time [s]')
ylabel('[A]')
title('AC currents')
legend('iL1', 'iL2', 'iL3')
grid on

aX2 = subplot(3,1,2);
plot(t, uDC);
xlabel('time [s]')
ylabel('[V]')
title('DC Voltage')
legend('uDC')
grid on

aX3 = subplot(3,1,3);
plot(t, iDC, 'r');
xlabel('time [s]')
ylabel('[A]')
title('DC Current')
legend('iDC')
grid on

linkaxes([aX3,aX2,aX1],'x');

%%% plot with 2 subplots
% aX1 = subplot(2,1,1);
% plot(t, [iL1 iL2 iL3])
% xlabel('time [s]')
% ylabel('[A]')
% title('AC currents')
% legend('iL1', 'iL2', 'iL3')
% grid on
% 
% aX2 = subplot(2,1,2);
% [hAx,hLine1,hLine2] = plotyy(t, uDC, t, iDC);
% title('DC')
% set(hLine1,'Color','blue')
% set(hLine2,'Color','red')
% set(hAx(2),'YColor','red')
% 
% % allign both axis with 0
% ylim1 = get(hAx(1),'Ylim');
% ylim2 = get(hAx(2),'Ylim');
% delta = ylim2(1) - ylim1(1)*((ylim2(2)-ylim2(1))/(ylim1(2)-ylim1(1)));
% set(hAx(2),'Ylim',[ylim2(1)-delta ylim2(2)-delta])
% 
% ylabel(hAx(1),'[V]') % left y-axis
% ylabel(hAx(2),'[A]') % right y-axis
% linkaxes([hAx(1), hAx(2)],'x');
% xlabel('time [s]')
% legend('uDC', 'iDC')
% grid on
% linkaxes([aX2,aX1],'x');

