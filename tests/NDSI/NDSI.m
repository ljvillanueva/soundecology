%Original Matlab code
% providede courtesy of S Gage on Oct 2013

%************************************************************
% Calculate the 'Average power spectral density (PSD) for each 1 kHz
% frequency band using Welch method.
%************************************************************
% nfft should be a power of 2
nfft = 1024;
Pxx = pwelch(signal,nfft,nfft/2,nfft,Fs);
Hpsd = dspdata.psd(Pxx,'Fs',Fs,'spectrumType','onesided'); % Create the object
%************************************************************
% Create the table matrix based on parameters set for acoustic spectrum
% analysis. Note that if a specific time is set, only one line of
% results will be generated instead of an entire table
%************************************************************

% Frequency bins work array
freqbins = zeros(1,10);

% Calculate the average power for frequency bins 1kHz to 11kHz
for k = 1:10
freqbins(k) = avgpower(Hpsd,[(k*1000) ((k+1)*1000)]);
end

% Vector normalize the PSD values between 0-1
freqbins = freqbins / norm(freqbins);
% Put normalized values into table, and set the first frequency bin
% to -1, indicating that the first 0-1 kHz is not computed.
Tab2(i,1) = -1;
for k = 1:10
Tab2(i,k+1) = freqbins(k);
end
% Calculate the average power for the entire signal from 1kHz to 11kHz.
% Sum the bins so that the total is on the same scale as the individual
% bins.
Tab2(i,12) = sum(freqbins);

% Get Anthrophony (1-2 bins)
Tab2(i,13) = Tab2(i,2);
% Get Biophony (3-11 bins) % Modified by SHG Oct 2012
Tab2(i,14) =
Tab2(i,3)+Tab2(i,4)+Tab2(i,5)+Tab2(i,6)+Tab2(i,7)+Tab2(i,8)+Tab2(i,9)+Tab2(i,10);T
ab2(i,11);
% Compute Normalized Difference Sound Index (NDSI B-A/B+A)
Tab2(i,15) = (Tab2(i,14)-Tab2(i,13))/(Tab2(i,14)+Tab2(i,13));
