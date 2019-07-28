function [data, truth] = helperCommsWaveformGeneration()

numFramesPerModType = 2;
sps = 8;                % Samples per symbol
spf = 1024;             % Samples per frame
fs = 100e6;             % Sampling rate

%% multipathChannel impairments

%AWGN
% snrVector = -20:5:20; % Range of Signal-to-
% Noise Ratios (SNRs) for training (dB)
SNR = 20;

% %CombinedmultipathChannel
% multipathChannel = helperModClassTestChannel(...
%   'SampleRate', fs, ...
%   'SNR', SNR, ...
%   'PathDelays', [0 1.8 3.4] / fs, ...
%   'AveragePathGains', [0 -2 -10], ...
%   'KFactor', 4, ...
%   'MaximumDopplerShift', 4, ...
%   'MaximumClockOffset', 5, ...
%   'CenterFrequency', 900e6);

multipathChannel = comm.RicianChannel(...
    'SampleRate', fs, ...
    'PathDelays', [0 1.8 3.4]/fs, ...
    'AveragePathGains', [0 -2 -10], ...
    'KFactor', 4, ...
    'MaximumDopplerShift', 4);

hFreqOffset = comm.PhaseFrequencyOffset(...
    'SampleRate',fs);

modulationTypes = categorical(["GFSK","CPFSK","2FSK","8FSK","BPSK","QPSK","8PSK","16QAM","64QAM","PAM4"]);

numModulationTypes = length(modulationTypes);

% Modified code by Rachana
idx = 1;
filterCoeffs = rcosdesign(0.35, 4, sps);
for iM = 1:numModulationTypes
     modType = modulationTypes(iM);
     rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
     switch modType
         case 'BPSK'
             for is = 1:numFramesPerModType
                 M = 2;
                 src = randi([0 M-1],spf/sps,1);
                 % bpskModulator BPSK modulator with pulse shaping
                 
                 % Modulate
                 syms = pskmod(src,2);
                 
                 % Pulse shape
                 y = filter(filterCoeffs, 1, upsample(syms,sps));
                                  
                 % Adjust SNR
                 wav = awgn(y,SNR);
                
                 % Add frequency offset
                 fc = randOverInterval(rangeFc);
                 hFreqOffset.FrequencyOffset = fc;
                 wav = hFreqOffset(wav); % Frequency shift
                 
                 rxSamples = multipathChannel(wav);
                 % Save signal
                 data{idx} = rxSamples;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
                 
                 release(hFreqOffset);
             end
         case 'QPSK'
             for is = 1:numFramesPerModType
                 M = 4;
                 src = randi([0 M-1],spf/sps,1);

                 % Modulate
                 syms = pskmod(src,4,pi/4);
                 % Pulse shape
                 y = filter(filterCoeffs, 1, upsample(syms,sps));
                 % Adjust SNR
                 wav = awgn(y,SNR);
                
                 % Add frequency offset
                 fc = randOverInterval(rangeFc);
                 hFreqOffset.FrequencyOffset = fc;
                 wav = hFreqOffset(wav); % Frequency shift
                 
                 rxSamples = multipathChannel(wav);
                 % Save signal
                 data{idx} = rxSamples;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
                 
                 release(hFreqOffset);
             end
         case '8PSK'
             for is = 1:numFramesPerModType
                 M = 8;
                 src = randi([0 M-1],spf/sps,1);

                 % Modulate
                 syms = pskmod(src,8);
                 % Pulse shape
                 y = filter(filterCoeffs, 1, upsample(syms,sps));
                % Adjust SNR
                 wav = awgn(y,SNR);
                
                 % Add frequency offset
                 fc = randOverInterval(rangeFc);
                 hFreqOffset.FrequencyOffset = fc;
                 wav = hFreqOffset(wav); % Frequency shift
                 
                 rxSamples = multipathChannel(wav);
                 % Save signal
                 data{idx} = rxSamples;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
                 
                 release(hFreqOffset);
             end
         case '2FSK'
             for is = 1:numFramesPerModType
                 M = 2;
                 freqsep = 4;  % Frequency separation (Hz)
                 symbolRate = 45;

                 src = randi([0 M-1],spf/sps,1);

                 % Modulate

                 fmod = comm.FSKModulator(...
                       'ModulationOrder', M, ...
                       'FrequencySeparation', freqsep, ...
                       'SymbolRate', symbolRate, ...
                       'SamplesPerSymbol', sps);
                 syms = fmod(src);
                 % Pulse shape
                 y = filter(filterCoeffs, 1, upsample(syms,sps));
                 % Adjust SNR
                 wav = awgn(y,SNR);
                
                 % Add frequency offset
                 fc = randOverInterval(rangeFc);
                 hFreqOffset.FrequencyOffset = fc;
                 wav = hFreqOffset(wav); % Frequency shift
                 
                 rxSamples = multipathChannel(wav);
                 % Save signal
                 data{idx} = rxSamples;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
                 release(fmod);
                 release(hFreqOffset);
             end
             clear fmod;
         case '8FSK'
             for is = 1:numFramesPerModType
                 M = 8;
                 src = randi([0 M-1],spf/sps,1);

                 % Modulate
                 syms = fskmod(src,M,freqsep,sps,fs);
                 % Pulse shape
                 y = filter(filterCoeffs, 1, upsample(syms,sps));
                 % Adjust SNR
                 wav = awgn(y,SNR);
                
                 % Add frequency offset
                 fc = randOverInterval(rangeFc);
                 hFreqOffset.FrequencyOffset = fc;
                 wav = hFreqOffset(wav); % Frequency shift
                 
                 rxSamples = multipathChannel(wav);
                 % Save signal
                 data{idx} = rxSamples;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
                 
                 release(hFreqOffset);
             end
             
         case 'GFSK'
             for is = 1:numFramesPerModType
                 M = 2;
                 src = randi([0 M-1],spf/sps,1);
                  
                 gmod = comm.CPMModulator(...
                    'ModulationOrder', M, ...
                    'FrequencyPulse', 'Gaussian', ...
                    'BandwidthTimeProduct', 0.35, ...
                    'ModulationIndex', 1, ...
                    'SamplesPerSymbol', sps);
                 meanM = mean(0:M-1);

                 % Modulate
                 y = gmod(2*(src-meanM));
                 
                 % Adjust SNR
                 wav = awgn(y,SNR);
                
                 % Add frequency offset
                 fc = randOverInterval(rangeFc);
                 hFreqOffset.FrequencyOffset = fc;
                 wav = hFreqOffset(wav); % Frequency shift
                 
                 rxSamples = multipathChannel(wav);
                 % Save signal
                 data{idx} = rxSamples;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
                 release(gmod);
                 release(hFreqOffset);
             end
             clear gmod;
         case '16QAM'
             for is = 1:numFramesPerModType
                 M = 16;
                 src = randi([0 M-1],spf/sps,1);

                 % Modulate and pulse shape
                 syms = qammod(src,16,'UnitAveragePower',true);
                 % Pulse shape
                 y = filter(filterCoeffs, 1, upsample(syms,sps));
                 % Adjust SNR
                 wav = awgn(y,SNR);
                
                 % Add frequency offset
                 fc = randOverInterval(rangeFc);
                 hFreqOffset.FrequencyOffset = fc;
                 wav = hFreqOffset(wav); % Frequency shift
                 
                 rxSamples = multipathChannel(wav);
                 % Save signal
                 data{idx} = rxSamples;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
                 
                 release(hFreqOffset);
             end
            
         case '64QAM'
             for is = 1:numFramesPerModType
                 M = 64;
                 src = randi([0 M-1],spf/sps,1);

                 % Modulate and pulse shape
                 syms = qammod(src,64,'UnitAveragePower',true);
                 % Pulse shape
                 y = filter(filterCoeffs, 1, upsample(syms,sps));
                 % Adjust SNR
                 wav = awgn(y,SNR);
                
                 % Add frequency offset
                 fc = randOverInterval(rangeFc);
                 hFreqOffset.FrequencyOffset = fc;
                 wav = hFreqOffset(wav); % Frequency shift
                 
                 rxSamples = multipathChannel(wav);
                 % Save signal
                 data{idx} = rxSamples;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
                 
                 release(hFreqOffset);
             end
         case 'CPFSK'
             for is = 1:numFramesPerModType
                 M = 2;
                 src = randi([0 M-1],spf/sps,1);

                   cpmod = comm.CPFSKModulator(...
                         'ModulationOrder', M, ...
                         'ModulationIndex', 0.5, ...
                         'SamplesPerSymbol', sps);
                   meanM = mean(0:M-1);

                 % Modulate
                 y = cpmod(2*(src-meanM));
                 % Adjust SNR
                 wav = awgn(y,SNR);
                
                 % Add frequency offset
                 fc = randOverInterval(rangeFc);
                 hFreqOffset.FrequencyOffset = fc;
                 wav = hFreqOffset(wav); % Frequency shift
                 
                 rxSamples = multipathChannel(wav);
                 % Save signal
                 data{idx} = rxSamples;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
                 release(cpmod);
                 release(hFreqOffset);
             end
             clear cpmod;
         case 'PAM4'
             for is = 1:numFramesPerModType
                 M = 4;
                 amp = 1 / sqrt(mean(abs(pammod(0:3, 4)).^2));
                 src = randi([0 M-1],spf/sps,1);

                 % Modulate
                 syms = amp * pammod(src,4);
                 % Pulse shape
                 y = filter(filterCoeffs, 1, upsample(syms,sps));
                 % Adjust SNR
                 wav = awgn(y,SNR);
                
                 % Add frequency offset
                 fc = randOverInterval(rangeFc);
                 hFreqOffset.FrequencyOffset = fc;
                 wav = hFreqOffset(wav); % Frequency shift
                 
                 rxSamples = multipathChannel(wav);
                 % Save signal
                 data{idx} = rxSamples;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
                 
                 release(hFreqOffset);
             end
     end
end

end

%% Subroutines
function val = randOverInterval(interval)
% Expect interval to be <1x2> with format [minVal maxVal]
val = (interval(2) - interval(1)).*rand + interval(1);
end

