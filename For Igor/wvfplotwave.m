function wvfplotwave(filename)
% plot waveforms of single .wvf file

[uL1,t] = wvfread(filename, 1, 1);
iL1 = wvfread(filename, 1, 2);
uL2 = wvfread(filename, 1, 3);
iL2 = wvfread(filename, 1, 4);
uL3 = wvfread(filename, 2, 1);
iL3 = wvfread(filename, 2, 2);
uDC = wvfread(filename, 2, 3);
iDC = wvfread(filename, 2, 4);

info = hdrread(filename);
datestring = [info.Group1.Trace1.Date{1,1}, ' ' info.Group1.Trace1.Time{1,1}];

% scaling correction
uL3 = 2*uL3;
uDC = 2*uDC;

% plots 
figure('Name', 'AC Voltages')
clf
plotbrowser('on')
plot(t,[uL1 uL2 uL3])
xlabel('[s]')
ylabel('[V]')
title(datestring)
legend('uL1', 'uL2', 'uL3')
grid on

figure('Name', 'AC Currents')
clf
plotbrowser('on')
plot(t,[iL1 iL2 iL3])
legend('iL1', 'iL2', 'iL3')
xlabel('[s]')
ylabel('[A]')
title(datestring);
grid on

figure('Name', 'DC')
clf
plotbrowser('on')


aX1 = subplot(2,1,1);
plot(t, uDC, 'b')
xlabel('[s]')
ylabel('[V]')
title(datestring)
legend('uDC')
grid on

aX2 = subplot(2,1,2);
plot(t, iDC, 'r')
xlabel('[s]')
ylabel('[A]')
legend('iDC')
grid on

linkaxes([aX2,aX1],'x');