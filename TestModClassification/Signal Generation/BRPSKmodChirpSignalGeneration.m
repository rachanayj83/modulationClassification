function [sentSignal, timeVec, RRCtransSignal, complexChirp, freqSymbol] = BRPSKmodChirpSignalGeneration(sampleTime, desDurTime, desSymFreq, rollOff, span, redPhase, startFreq, stopFreq, chirpPhase, chirpAmp)
%
% [sentSignal, timeVec, RRCtransSignal, complexChirp, freqSymbol] = 
% BRPSKmodChirpSignalGeneration(sampleTime, desDurTime, desSymFreq, 
% rollOff, span, redPhase, startFreq, stopFreq, chirpPhase, chirpAmp)
%
% Generates a BRPSK modulated chirp signal, take into account that if the 
% reduced phase (redPhase) is equal to the 90 degrees, or pi/2 radians, the
% sent signal is the same as a BPSK modulated chirp signal.
%
% INPUTS:
% sampleTime    - Sample time
% desDurTime    - Desired duration of the chirp signal 
% desSymFreq    - Desired symbol frequency (as the data is binary, the 
%               symbol and the bit frequency are the same)
% rollOff       - Roll off factor of the RRC filter
% span          - Span for the RRC filter
% redPhase      - Reduced phase of the BRPSK modulation
% startFreq     - Start frequency of the chirp signal
% stopFreq      - Stop frequenncy of the chirp signal
% chirpPhase    - Initial phase of the chirp 
% chirpAmp      - Amplitude of the chirp
% 
% OUTPUTS:
% sentSignal        - Sent signal
% timeVec           - Time vector 
% RRCtransSignal    - Output signal of the RRC filter
% complexChirp      - Complex envelope chirp signal
% freqSymbol        - Actual symbol frequency
%

%% Simulation data
[numSymbols, numSampPerSym, freqSymbol, ~, ~] = simDataCalc(desDurTime, desSymFreq, sampleTime);

%% RRC filter characteristics
hTxFilter = comm.RaisedCosineTransmitFilter('RolloffFactor',rollOff,...
'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',numSampPerSym);

coefficients = coeffs(hTxFilter);
hTxFilter.Gain = 1/sum(coefficients.Numerator);

%% n Guard Bits
B = 1;
nG = ceil(span*B);

%% LFM chirp signal
chirpTime = ((numSymbols + nG)*numSampPerSym - 1)*sampleTime;

complexChirp = complexChirpSignalGenerator(startFreq, stopFreq, chirpAmp, chirpPhase, chirpTime, sampleTime);

%% Data generation
seqInf = binDataGeneration(numSymbols);
G = binDataGeneration(nG);
binArray = [seqInf; G];

%% BRPSK modulation
modSymbols = BRPSKsymbols(binArray, redPhase);

%% RRC transmit filter
reset(hTxFilter);
RRCtransSignal = step(hTxFilter, modSymbols);

sentSignal = real(complexChirp.*RRCtransSignal);

powerSignal = (sum(abs(sentSignal).^2)/length(sentSignal));
sentSignal = sentSignal/sqrt(powerSignal);

timeVec = (0:(length(sentSignal) - 1))*sampleTime;
end