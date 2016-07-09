%% mean values calculation

i_Eff = [2]; % column index of 'X' and 'header' for efficiency
i_Pow = [3 4 141]; % column indices of 'X' and 'header' for PAC and PDC
i_QcosphiPF = [5 6 9 10]; % column indices of 'X' and 'header' for Q1, Q1D, cosphi and PF
i_UL = [11 12 13]; % column indices of 'X' and 'header' for rms values of AC voltages
i_IL = [14 15 16]; % column indices of 'X' and 'header' for rms values of AC currents
i_UIDC = [17 18]; % column indices of 'X' and 'header' for DC voltage and current
i_THDuharmL1 = [19 25:43]; % column indices of 'X' and 'header' for THDu and voltage harmonics up to order 20
i_THDiharmL1 = [22 44:62]; % column indices of 'X' and 'header' for THDi and current harmonics up to order 20
i_THDuharmL2 = [20 63:81]; % column indices of 'X' and 'header' for THDu and voltage harmonics up to order 20
i_THDiharmL2 = [23 82:100]; % column indices of 'X' and 'header' for THDi and current harmonics up to order 20
i_THDuharmL3 = [21 101:119]; % column indices of 'X' and 'header' for THDu and voltage harmonics up to order 20
i_THDiharmL3 = [24 120:138]; % column indices of 'X' and 'header' for THDi and current harmonics up to order 20
i_ripple = [139 140];  % column indices of 'X' and 'header' for voltage ripple and current ripple

lim = [
    14 41;  % lim for 10 kW (row indeces of Xread)
    43 57;  % lim for 7.5kW (row indeces of Xread)
    60 78;  % lim for 5 kW  (row indeces of Xread)
    80 117; % ...
    121 141;
    143 165;
    167 183;
    185 218;
    ];

% choose which values you want to average!
indeces = [i_Eff i_Pow i_THDuharmL3 i_THDiharmL3 i_ripple];

Av = zeros(size(lim,1),length(indeces));
for i=1:size(lim,1)
    Av(i,:) = mean(cell2mat(Xread(lim(i,1):lim(i,2),indeces)));
end
Av = [header(indeces); num2cell(Av)]; 
% disp(Av);