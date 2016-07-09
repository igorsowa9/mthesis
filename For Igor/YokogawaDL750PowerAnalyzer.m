%% Yokogawa DL750 Power Analyzer
% For stored data by "action on trigger" in binary format (.wvf and .hdr)
% The channels should have these assignments:
%
% Ch1: uL1
% Ch2: iL1
% Ch3: uL2
% Ch4: iL2
% Ch5: uL3
% Ch6: iL3
% Ch7: uDC
% Ch8: iDC
% 
% If this is not the case, one has to change the assignments in the
% wvfcalcvalues.m function and maybe change the plots commands in section 2.
%
% It is also possible to change the calculated values in wvfcalcvalues.m
% and adapt them in this script.
%
% Needs the functions wvfread.m, wvfreadb.m, hdrread.m, wvfcalcvalues.m,
% wvfplotwave.m, dftgeneral.m and harm50.m in the same folder.

%% Section1: Calculate desired values out of set of binary files (.wvf and .hdr)
% The scirpt automatically saves fnames, header, X and Xread once
% executed. You can then just load the
% matlab variables, skip this section and execute the following sections with the
% loaded data.

[fnames, pname] = uigetfile('*.wvf', 'Select WVF files', 'Multiselect', 'on');
header = wvfcalcvalues('header');
numFiles = numel(fnames);
numHeader = numel(header);

X = zeros(numFiles, numHeader);    
for i = 1:numFiles
    X(i,:) = wvfcalcvalues([pname fnames{i}]);
end

Xread = [cellstr(datestr(X(:,1), 'dd.mm.yy HH:MM:SS')), num2cell(X(:, 2:end))];
Xread = [
    header;
    Xread
    ];
save([pname fnames{1} '-' fnames{end} '.mat'], 'X', 'Xread', 'pname', 'fnames', 'header');

%% Section 2a: Plots

% efficiencies
i_Eff = [2]; % column index of 'X' and 'header' for efficiency
figure(1)
set(gcf, 'Name', 'Efficiency')
clf
plotbrowser('on')
plot(X(:,1), X(:,i_Eff))
hold on
grid on
legend(header{i_Eff}) 
datetick
xlabel('Date and Time')
ylabel('Efficiency [%]')

% powers
i_Pow = [3 4 141]; % column indices of 'X' and 'header' for PAC, PDC and Plosses
figure(2)
set(gcf, 'Name', 'Powers')
clf
plotbrowser('on')
plot(X(:,1), X(:,i_Pow))
grid on
legend(header{i_Pow})
datetick
xlabel('Date and Time')
ylabel('Power [W]')

% Q, cosphi and PF
i_QcosphiPF = [5 6 9 10]; % column indices of 'X' and 'header' for Q1, Q1D, cosphi and PF
figure(3)
set(gcf, 'Name', 'Reactive powers and power factors')
clf
plotbrowser('on')

aX1 = subplot(2,1,1);
plot(X(:,1), X(:,i_QcosphiPF(1:2)))
grid on
legend(header{i_QcosphiPF(1:2)})
xlabel('Date and Time')
ylabel('Reactive Power [var]')
datetick

aX2 = subplot(2,1,2);
plot(X(:,1), X(:,i_QcosphiPF(3:4)))
grid on
legend(header{i_QcosphiPF(3:4)})
xlabel('Date and Time')
ylabel('[1]')
datetick

linkaxes([aX2,aX1],'x');

% grid voltages
i_UL = [11 12 13]; % column indices of 'X' and 'header' for rms values of AC voltages
figure(4)
set(gcf, 'Name', 'Grid Voltages')
clf
plotbrowser('on')
plot(X(:,1), X(:,i_UL))
grid on
legend(header{i_UL})
datetick
xlabel('Date and Time')
ylabel('Voltage [V]')

% grid currents
i_IL = [14 15 16]; % column indices of 'X' and 'header' for rms values of AC currents
figure(5)
set(gcf, 'Name', 'Grid Currents')
clf
plotbrowser('on')
plot(X(:,1), X(:,i_IL))
grid on
legend(header{i_IL})
datetick
xlabel('Date and Time')
ylabel('Current [A]')

% DC voltage and current
i_UIDC = [17 18]; % column indices of 'X' and 'header' for DC voltage and current
figure(6)
set(gcf, 'Name', 'DC')
clf
plotbrowser('on')

ax1 = subplot(2,1,1);
plot(X(:,1), X(:,i_UIDC(1)))
hold on
grid on
legend(header{i_UIDC(1)})
xlabel('Date and Time')
ylabel('Voltage [V]')
datetick

ax2 = subplot(2,1,2);
plot(X(:,1), X(:,i_UIDC(2)), 'r')
hold on
grid on
legend(header{i_UIDC(2)})
xlabel('Date and Time')
ylabel('Current [A]')
datetick

linkaxes([ax2,ax1],'x');

% THD and odd harmonics trend U
i_THDuoddharmL1 = [19 26:2:42]; % column indices of 'X' and 'header' for THDu and odd voltage harmonics L1
i_THDuoddharmL2 = [20 64:2:80]; % column indices of 'X' and 'header' for THDu and odd voltage harmonics L2
i_THDuoddharmL3 = [21 102:2:118]; % column indices of 'X' and 'header' for THDu and odd voltage harmonics L3
figure(7)
set(gcf, 'Name', 'Voltage Harmonics L1 Trend')
clf
plotbrowser('on')
cc = hsv(length(i_THDuoddharmL1));
hold on
for i=1:length(i_THDuoddharmL1)
    plot(X(:,1), X(:,i_THDuoddharmL1(i)), 'Color', cc(i,:));
end
grid on
legend(header{i_THDuoddharmL1})
datetick
xlabel('Date and Time')
ylabel('[%]')


% THD and odd harmonics trend I
i_THDioddharmL1 = [22 45:2:61]; % column indices of 'X' and 'header' for THDi and odd current harmonics L1
i_THDioddharmL2 = [23 83:2:99]; % column indices of 'X' and 'header' for THDi and odd current harmonics L2
i_THDioddharmL3 = [24 121:2:137]; % column indices of 'X' and 'header' for THDi and odd current harmonics L3

figure(8)
set(gcf, 'Name', 'Current Harmonics L1 Trend')
clf
plotbrowser('on')
cc = hsv(length(i_THDioddharmL1));
hold on
for i=1:length(i_THDioddharmL1)
    plot(X(:,1), X(:,i_THDioddharmL1(i)), 'Color', cc(i,:));
end
grid on
legend(header{i_THDioddharmL1})
datetick
xlabel('Date and Time')
ylabel('[%]')

% ripple of DC voltage and current
i_ripple = [139 140]; % column indices of 'X' and 'header' for DC voltage and current
figure(9)
set(gcf, 'Name', 'ripple')
clf
plotbrowser('on')

ax1 = subplot(2,1,1);
plot(X(:,1), X(:,i_ripple(1)))
hold on
grid on
legend(header{i_ripple(1)})
xlabel('Date and Time')
ylabel('Voltage Ripple [V]')
datetick

ax2 = subplot(2,1,2);
plot(X(:,1), X(:,i_ripple(2)), 'r')
hold on
grid on
legend(header{i_ripple(2)})
xlabel('Date and Time')
ylabel('Current Ripple [A]')
datetick

linkaxes([ax2,ax1],'x');

%% Section 2b: Harmonics Boxplots

% assigning the column indices
i_THDuharmL1 = [19 25:43]; % column indices of 'X' and 'header' for THDu and voltage harmonics up to order 20
i_THDiharmL1 = [22 44:62]; % column indices of 'X' and 'header' for THDi and current harmonics up to order 20
i_THDuharmL2 = [20 63:81]; % column indices of 'X' and 'header' for THDu and voltage harmonics up to order 20
i_THDiharmL2 = [23 82:100]; % column indices of 'X' and 'header' for THDi and current harmonics up to order 20
i_THDuharmL3 = [21 101:119]; % column indices of 'X' and 'header' for THDu and voltage harmonics up to order 20
i_THDiharmL3 = [24 120:138]; % column indices of 'X' and 'header' for THDi and current harmonics up to order 20


% getting start and end datetime and line
prompt = {'start datetime (dd.mm.yyyy HH:MM:SS):','end datetime (dd.mm.yyyy HH:MM:SS):', 'Line'};
dlg_title = 'Harmonics time frame';
num_lines = 1;
lim = xlim;
defaultans = {datestr(lim(1),'dd.mm.yyyy HH:MM:SS'),datestr(lim(2),'dd.mm.yyyy HH:MM:SS'), '1'};
input = inputdlg(prompt,dlg_title,num_lines,defaultans);
datetimelim = [datenum(input{1}, 'dd.mm.yyyy HH:MM:SS') datenum(input{2}, 'dd.mm.yyyy HH:MM:SS')];
line = input{3};
i_time = find(datetimelim(1)<X(:,1) & X(:,1)<datetimelim(2));

if strcmp(line, '1')
    i_U = i_THDuharmL1;
    i_I = i_THDiharmL1;
else if strcmp(line, '2')
    i_U = i_THDuharmL2;
    i_I = i_THDiharmL2;
else if strcmp(line, '3')
    i_U = i_THDuharmL3;
    i_I = i_THDiharmL3;
    end
    end
end

% plot
figure()
set(gcf, 'Name', ['Boxplots Harmonics L', line])
clf
plotbrowser('on')

AX1 = subplot(2,1,1);
boxplot(X(i_time, i_U),'Labels', header(i_U),'labelorientation', 'inline')
grid on
title(['Voltage Harmonics from ', input{1}, ' to ', input{2}]);
ylabel('[%]')

AX2 = subplot(2,1,2);
boxplot(X(i_time, i_I),'Labels', header(i_I),'labelorientation', 'inline')
grid on
title(['Current Harmonics from ', input{1}, ' to ', input{2}]);
ylabel('[%]')

%% Section 3: Possibility to plot waveform of current cursor position
% Use the cursor in a desired plot, choose a point, right click and save to
% workspace variable named 'cursor_info'. Then execute this section.
wvfplotwave([pname fnames{cursor_info.DataIndex}]);
