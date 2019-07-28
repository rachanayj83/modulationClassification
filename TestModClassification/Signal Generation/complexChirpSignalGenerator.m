function complexChirp = complexChirpSignalGenerator(startFreq, stopFreq, complexAmp, startPhase, chirpDuration, sampleTime)
%
% complexChirp = complexChirpSignalGenerator(startFreq, stopFreq, 
% complexAmp, startPhase, chirpDuration, sampleTime)
% Generates a complex up-chirp signal
%
% INPUTS:
% startFreq     - start frequency of the chirp
% stopFreq      - stop frequency of the chirp
% complexAmp    - complex amplitude of the chirp signal
% startPhase    - start phase of the complex chirp signal
% chirpDuration - duration of the chirp (in seconds)
% sampleTime    - duration of one sample of the chirp signal
%
% OUTPUTS:
% complexChirp  - complex up chirp signal with the desired characteristics
%

k = (stopFreq - startFreq)/chirpDuration;

tChirp = transpose(0:sampleTime:chirpDuration);

complexChirp = complexAmp*exp(1i*(2*pi*(startFreq + (k/2)*tChirp).*tChirp + startPhase));