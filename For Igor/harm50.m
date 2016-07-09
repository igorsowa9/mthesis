function [HARM50] = harm50(FFT)
% calculates the harmonics from FFT matrix (which can be obtained via the
%   dftgeneral function) in accordance with IEEE Std 519-2014.
%   The frequency stepsize of the FFT matrix needs to be 5 Hz. This will be
%   accomplished by a measurement over 200 ms (10 cycles) and/or by
%   choosing the right number of points from the waveformvector.
%   The HARM50 matrix has the following form:
%
%   HARM50 = [  frequency of harmonic (Hz);
%               rms magnitude (unit of signal);
%               proportion to fundamental (%);
%               phase (º);                     ]
%
%   Moreover, the column index indicates the order of the harmonic.
%   The phases are calculated for cosine and are shifted so that the phase
%   of the fundamental is 0.

FFT(1,:) = int64(FFT(1,:)); % to avoid number format that could lead to problems.

if (FFT(1,2)-FFT(1,1)) ~= 5
    errordlg(['Not possible to determine harmonics because the frequency stepsize in the FFT matrix ',...
    'is not 5 Hz. The measurement has to be done over 10 cycles (200ms). You could also adjust ',...
    'the number of points you take from your waveformvector.'],'Wrong frequency stepsize')
end

ind = find(mod(FFT(1,1:end),50)==0);
num_harm = length(ind)-1;
fund_rms = sqrt(sumsqr(FFT(2,ind(2)-1:ind(2)+1))/2);
fund_phi = FFT(3,ind(2));

HARM50 = zeros(4,num_harm);

for i = 1:num_harm
    harm_freq = FFT(1, ind(i+1));
    harm_rms = sqrt(sumsqr(FFT(2,ind(i+1)-1:ind(i+1)+1))/2);
    harm_phase = mod(FFT(3, ind(i+1))-i*fund_phi, 360); % Phase Shift, so that fundamental phase is zero
    if harm_phase > 180; harm_phase = harm_phase - 360; end
    HARM50(:,i)=[
        harm_freq;
        harm_rms;
        harm_rms/fund_rms*100;
        harm_phase;
        ];
end

end