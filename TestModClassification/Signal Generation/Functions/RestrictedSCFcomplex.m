function [resSCF, SCF, ACF1, freq, alphaFreq, limitMatrix] = RestrictedSCFcomplex(signal, alphaResol, alphaLimit, freqSCFresol, sampleTime, restrictFactor)
%
% [resSCF, SCF, ACF1, freq, alphaFreq, limitMatrix] = RestrictedSCFcomplex
% (signal, alphaResol, alphaLimit, freqSCFresol, sampleTime, 
% restrictFactor)
% 
% Computes the restricted complex Spectral Correlation Function of the
% input signal with the alpha frequency resolution, alphaResol, the SCF
% frequency resolution, freqSCFresol.
%
% It must be said that the freqSCFresol, the SCF frequency resolution, 
% might change a little for assorting the computing operations.
%
% INPUTS:
% signal                - Input signal
% alphaResol            - Resolution of the alpha frequencies
% alphaLimit            - Up-limit of the alpha frequencies
% freqSCFresol          - Resolution of the SCF frequencies
% sampleTime            - Sample time
% restrictFactor        - Restrictive factor in decimals
% 
% OUTPUTS:
% resSCF                - Restricted Spectral Correlation Function of the
%                       signal
% SCF                   - Non-restricted Spectral Correlation Function
% freq                  - Vector of frequencies
% alphaFreq             - Vector of the alpha frequencies
% limitMatrix           - Matrix with the limits that have been used for
%                       computing the restricted SCF
% 

numSamples = length(signal);
signalTime = (numSamples - 1)*sampleTime;

alphaFreq = transpose((-1*alphaLimit):alphaResol:alphaLimit);

tauResol = sampleTime;
freqSCFsample = 1/tauResol;
numFreqSCFsamples = round(freqSCFsample/freqSCFresol);
sampleTauLimit = round(numFreqSCFsamples/2);

tau = transpose((-1*sampleTauLimit*tauResol):tauResol:(sampleTauLimit*tauResol));
tauSamples = round(tau/sampleTime);

freqResol = freqSCFsample/length(tau); % Actual frequency resolution the function is working with
freqLimit = floor(length(tau)/2)*freqResol;
freq = transpose(-freqLimit:freqResol:freqLimit);

timeVec = transpose(0:sampleTime:signalTime);

numZerosTime1 = numSamples - 1;
funcTime1 = [zeros(numZerosTime1, 1); signal; zeros(numZerosTime1, 1)];

startFuncTime = numSamples;
stopFuncTime = 2*numSamples - 1;

ACF1 = zeros(length(alphaFreq), length(tau));

for iACF1 = 1:length(alphaFreq)    
    for iACF2 = 1:length(tau)
        numZeros2 = round(numSamples - 1 + tauSamples(iACF2));
        numZeros3 = round(numSamples - 1 - tauSamples(iACF2));
        
        funcTime2 = [zeros(numZeros2, 1); signal; zeros(numZeros3, 1)];
        
        multiplication = funcTime1(startFuncTime:stopFuncTime).*conj(funcTime2(startFuncTime:stopFuncTime));
        ACF1(iACF1, iACF2) = (1/signalTime)*(sum(multiplication.*exp(-1i*2*pi*alphaFreq(iACF1)*timeVec))*exp(1i*pi*alphaFreq(iACF1)*tau(iACF2)))*sampleTime;%**CAMBIAR ESTO****sampleTime))

    end
end

%% Restricted ACF
tauOmega = tau(round((length(tau)/2)*restrictFactor));
tauGood = (abs(tau) < abs(tauOmega));
multiCol = ones(length(alphaFreq), 1);

limitMatrix = multiCol*transpose(tauGood);
ACF2 = ACF1.*limitMatrix;

%% Restricted SCF
SCF = zeros(length(alphaFreq), length(tau));
resSCF = zeros(length(alphaFreq), length(tau));

for iSCF1 = 1:length(alphaFreq)
    SCF(iSCF1, :) = fftshift(fft(ACF1(iSCF1, :)));
    resSCF(iSCF1, :) = fftshift(fft(ACF2(iSCF1, :)));
end

end