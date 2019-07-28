function [numSymbols, sampPerSymbol, freqSymbol, numSamples, finSignalTime] = simDataCalc(desiredSignalTime, desiredSymbolFreq, sampleTime)
%
% [numSymbols, sampPerSymbol, freqSymbol, numSamples, finSignalTime] = 
% simDataCalc(desiredSignalTime, desiredSymbolFreq, sampleTime)
%
% Computes the data for the simulation.
% 
% INPUTS:
% desiredSignalTime     - desired time of the final signal
% desiredSymbolFreq     - desired symbol frequency
% sampleTime            - time at which the signal is being sampled
%
% OUTPUTS:
% numSymbols        - number of symbols
% sampPerSymbol     - number of samples per symbol
% freqSymbol        - final symbol frequency
% numSamples        - number of samples
% finSignalTime     - final signal time
%

symTime = 1/desiredSymbolFreq;
numSymbols = round((desiredSignalTime + sampleTime)/symTime);

sampPerSymbol = round(symTime/sampleTime);

numSamples = numSymbols*sampPerSymbol;

finSymTime = sampPerSymbol*sampleTime;
freqSymbol = 1/finSymTime;

finSignalTime = (numSamples - 1)*sampleTime;
end