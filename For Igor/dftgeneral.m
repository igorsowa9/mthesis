function [FFTmatrix] = dftgeneral(waveformvector, Fs)
% executes useful DFT from a waveformvector with known sample frequency Fs
%   FFTmatrix has the following structure:
%       [   frequency(Hz); 
%           amplitude(unit of signal);
%           angle(º)     ]

x=waveformvector;
nfft = length(x);

X = fft(x);
X = X(1:nfft/2);  % FFT is symmetric, throw away the other half

f = (0:nfft/2-1)*Fs/nfft;  % Frequency vector
mx = abs(X*2/nfft); % Magnitude (Amplitude)
mx(1) = mx(1)/2;    % Offset (f=0) divided by 2
ax = angle(X)*180/pi; % Angle in º

% create desired FFTmatrix
FFTmatrix = zeros(3, nfft/2);
FFTmatrix(1, :) = f;
FFTmatrix(2, :) = mx;
FFTmatrix(3, :) = ax;

end