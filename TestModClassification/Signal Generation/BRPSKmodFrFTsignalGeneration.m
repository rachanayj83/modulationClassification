function [sentSignal, timeVec, RRCtransSignal, freqSymbol] = BRPSKmodFrFTsignalGeneration(sampleTime, desDurTime, desSymFreq, rollOff, span, redPhase, order)
%
% [sentSignal, timeVec, RRCtransSignal, freqSymbol] = 
% BRPSKmodFrFTsignalGeneration(sampleTime, desDurTime, desSymFreq, rollOff,
% span, redPhase, order)
%
% Generates a BRPSK modulated FrFT signal, take into account that if the 
% reduced phase (redPhase) is equal to the 90 degrees, or pi/2 radians, the
% sent signal is the same as a BPSK modulated FrFT signal.
%
% INPUTS:
% sampleTime    - Sample time
% desDurTime    - Desired duration of the chirp signal 
% desSymFreq    - Desired symbol frequency (as the data is binary, the 
%               symbol and the bit frequency are the same)
% rollOff       - Roll off factor of the RRC filter
% span          - Span for the RRC filter
% redPhase      - Reduced phase of the BRPSK modulation
% order         - Order of the FrFT
% 
% OUTPUTS:
% sentSignal        - Sent signal
% timeVec           - Time vector 
% RRCtransSignal    - Output signal of the RRC filter
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

%% Data generation
seqInf = binDataGeneration(numSymbols);
G = binDataGeneration(nG);
binArray = [seqInf; G];

%% BRPSK modulation
modSymbols = BRPSKsymbols(binArray, redPhase);

%% RRC transmit filter
reset(hTxFilter);
RRCtransSignalPre = step(hTxFilter, modSymbols);

numSamplesPre = length(RRCtransSignalPre);
RRCtransSignal = [RRCtransSignalPre; zeros(mod(numSamplesPre, 2) + 1, 1) + 1i*zeros(mod(numSamplesPre, 2) + 1, 1)];

%% FrFT
frftSignal = fracft(RRCtransSignal, order);

powerSignal = (sum(abs(frftSignal).^2)/length(frftSignal));
sentSignal = frftSignal/sqrt(powerSignal);

timeVec = (0:(length(sentSignal) - 1))*sampleTime;
end