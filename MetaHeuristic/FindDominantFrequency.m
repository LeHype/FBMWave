function [DF] = FindDominantFrequency(Signal,Fs)
%FINDDOMINANTFREQUENCY Returns the dominant frequency of signal
y = Signal;
NFFT = length(y);
Y = fft(y,NFFT);
F = ((0:1/NFFT:1-1/NFFT)*Fs).';
magnitudeY = abs(Y);
F = F(1:floor(end/2));
magnitudeY = magnitudeY(1:floor(end/2));

[a idx] = max(magnitudeY);
DF = 1/F(idx);
end

