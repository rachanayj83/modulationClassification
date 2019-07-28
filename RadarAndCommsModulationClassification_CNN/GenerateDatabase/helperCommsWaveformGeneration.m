function [data, truth] = helperCommsWaveformGeneration(frameStore)
%% User defined parameters
numFramesPerModType = 10;
sps = 8;                % Samples per symbol
spf = 512;             % Samples per frame
fs = 100e3;             % Sampling rate

%% Initialisation
modulationTypes = categorical(["LFM","Rect","Barker","Frank","P1","P2","P3","P4","Zadoff-Chu",...
    "BPSK","QPSK","8PSK","PAM4",...
    "16QAM","64QAM","GFSK","CPFSK","2FSK","8FSK"]);%

%  "B-FM", "DSB-AM", "SSB-AM"]);
%% Channel impairments
%AWGN
SNR = 10;
%snrVector = -20:5:20; % Range of Signal-to-Noise Ratios (SNRs) for training (dB)
%CombinedChannel
channel = helperModClassTestChannel(...
  'SampleRate', fs, ...
  'SNR', SNR, ...
  'PathDelays', [0 1.8 3.4] / fs, ...
  'AveragePathGains', [0 -2 -10], ...
  'KFactor', 4, ...
  'MaximumDopplerShift', 4, ...
  'MaximumClockOffset', 5, ...
  'CenterFrequency', 900e6);% Digital modulation types use 900 MHz center frequency
   transDelay = 50;
   
   multipathChannel = comm.RicianChannel(...
    'SampleRate', fs, ...
    'PathDelays', [0 1.8 3.4]/fs, ...
    'AveragePathGains', [0 -2 -10], ...
    'KFactor', 4, ...
    'MaximumDopplerShift', 4);

hFreqOffset = comm.PhaseFrequencyOffset(...
    'SampleRate',fs);
   
filterCoeffs = rcosdesign(0.35, 4, sps);
numModulationTypes = length(modulationTypes);

rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
fc = randOverInterval(rangeFc);

idx = 1;

for iM = 1:numModulationTypes
    modType = modulationTypes(iM);
     switch modType
         case 'BPSK'
             for is = 1:numFramesPerModType
                 M = 2;
                 src = randi([0 M-1],spf/sps,1);
                 % bpskModulator BPSK modulator with pulse shaping
                 % Y=bpskModulator(X,SPS) BPSK modulates input, X, and returns a root-raised
                 % cosine pulse shaped signal, Y. X must be a column vector of values in
                 % the set of [0 1]. Root-raised cosine filter has a roll-off factor of
                 % 0.35 and spans four symbols. The output signal, Y, has unit power.
                  if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                  end
                 % Modulate
                 syms = pskmod(src,2);
                
                 % Pass through independent channels
                 y = channel(syms);
                 
                  % Pulse shape
                 wav = filter(filterCoeffs, 1, upsample(y,sps));
                 
                 % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, spf, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modulationTypes(modType));
                 % Save signal
                 data{idx} = wav;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
             end
         case 'QPSK'
             for is = 1:numFramesPerModType
                 M = 4;
                 src = randi([0 M-1],spf/sps,1);
                 if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                  end
                 % Modulate
                 syms = pskmod(src,4,pi/4);
                 % Pass through independent channels
                 y = channel(syms);
                 
                  % Pulse shape
                 wav = filter(filterCoeffs, 1, upsample(y,sps));
                 
                  % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, spf, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modulationTypes(modType));
                 % Save signal
                 data{idx} = wav;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
             end
         case '8PSK'
             for is = 1:numFramesPerModType
                 M = 8;
                 src = randi([0 M-1],spf/sps,1);
                 if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                 end
                 % Modulate
                 syms = pskmod(src,8);
                 % Pass through independent channels
                 y = channel(syms);
                 
                  % Pulse shape
                 wav = filter(filterCoeffs, 1, upsample(y,sps));
                 
                  % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, spf, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modulationTypes(modType));
                 % Save signal
                 data{idx} = wav;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
             end
         case '2FSK'
             for is = 1:numFramesPerModType
                 M = 2;
                 freqsep = 8;  % Frequency separation (Hz)
                 symbolRate = 45;
%                  fsamp = sps*symbolRate;
                 src = randi([0 M-1],spf/sps,1);
                 if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                  end
                 % Modulate
%                   syms = fskmod(src,M,freqsep,sps,fs);
                 fmod = comm.FSKModulator(...
                       'ModulationOrder', M, ...
                       'FrequencySeparation', freqsep, ...
                       'SymbolRate', symbolRate, ...
                       'SamplesPerSymbol', sps);
                 syms = fmod(src);
                 
                 % Pulse shape
                 wav = filter(filterCoeffs, 1, upsample(syms,sps));
                 
                 % Adjust SNR
                wav = awgn(wav,SNR);
                
                % Add frequency offset
                hFreqOffset.FrequencyOffset = fc;
                wav = hFreqOffset(wav); % Frequency shift
                
                % Add multipath offset
                wav = multipathChannel(wav);

                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, spf, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modulationTypes(modType));
                 % Save signal
                 data{idx} = wav;
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
                 if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                 end
                 % Modulate
                 syms = fskmod(src,M,freqsep,sps,fs);
                  
                  % Pulse shape
                  wav = filter(filterCoeffs, 1, upsample(syms,sps));
                  
                   % Adjust SNR
                wav = awgn(wav,SNR);
                
                % Add frequency offset
                hFreqOffset.FrequencyOffset = fc;
                wav = hFreqOffset(wav); % Frequency shift
                
                % Add multipath offset
                wav = multipathChannel(wav);
 % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, spf, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modulationTypes(modType));
                 % Save signal
%                  rxSamples = awgn(y,SNR);
                 data{idx} = wav;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
                 release(hFreqOffset);
             end
         case 'GFSK'
             for is = 1:numFramesPerModType
                 M = 2;
                 src = randi([0 1],spf/sps,1);
                 %   Y=gfskModulator(X,SPS) GFSK modulates input, X, and returns signal, Y. 
                 %   X must be a column vector of values in the set of [0 1]. BT product is
                 %   0.35 and modulation index is 1. The output signal, Y, has unit power.

%                  if isempty(mod)
                 gmod = comm.CPMModulator(...
                    'ModulationOrder', M, ...
                    'FrequencyPulse', 'Gaussian', ...
                    'BandwidthTimeProduct', 0.35, ...
                    'ModulationIndex', 1, ...
                    'SamplesPerSymbol', sps);
                 meanM = mean(0:M-1);
%                  end
                 % Modulate
                 syms = gmod(2*(src-meanM));
                 
                 % Pass through independent channels
%                  wav = channel(syms);
                 
                 % Pulse shape
                 wav = filter(filterCoeffs, 1, upsample(syms,sps));
                 
                   % Adjust SNR
                wav = awgn(wav,SNR);
                
                % Add frequency offset
                hFreqOffset.FrequencyOffset = fc;
                wav = hFreqOffset(wav); % Frequency shift
                
                % Add multipath offset
                wav = multipathChannel(wav);
                 
                  % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, spf, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modulationTypes(modType));
                 % Save signal
                 data{idx} = wav;
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
                 if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                 end
                 % Modulate and pulse shape
                 syms = qammod(src,16,'UnitAveragePower',true);
                  % Pass through independent channels
                 y = channel(syms);
                 
                  % Pulse shape
                 wav = filter(filterCoeffs, 1, upsample(y,sps));
                 
                  % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, spf, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modulationTypes(modType));
                 % Save signal
                 data{idx} = wav;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
             end
         case '64QAM'
             for is = 1:numFramesPerModType
                 M = 64;
                 src = randi([0 M-1],spf/sps,1);
                 if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                 end
                 % Modulate and pulse shape
                 syms = qammod(src,64,'UnitAveragePower',true);
                 % Pass through independent channels
                 y = channel(syms);
                 
                  % Pulse shape
                 wav = filter(filterCoeffs, 1, upsample(y,sps));
                 
                  % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, spf, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modulationTypes(modType));
                 % Save signal
                 data{idx} = wav;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
             end
         case 'CPFSK'
             for is = 1:numFramesPerModType
                 M = 2;
%                   fc = randOverInterval(rangeFc);
                 src = randi([0 1],spf/sps,1);
%                  if isempty(mod)
                   cpmod = comm.CPFSKModulator(...
                         'ModulationOrder', M, ...
                         'ModulationIndex', 0.5, ...
                         'SamplesPerSymbol', sps);
                   meanM = mean(0:M-1);
%                  end
                 % Modulate
                 syms = cpmod(2*(src-meanM));
                 
                  wav = filter(filterCoeffs, 1, upsample(syms,sps));

                  % Adjust SNR
                  wav = awgn(wav,SNR);
                
                % Add frequency offset
                hFreqOffset.FrequencyOffset = fc;
                wav = hFreqOffset(wav); % Frequency shift
                
%                 Add multipath offset
                wav = multipathChannel(wav);
                  % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, spf, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modulationTypes(modType));
                 
                 % Save signal
                 data{idx} = wav;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
                 release(hFreqOffset);
                 release(cpmod);
             end
             clear cpmod;
         case 'PAM4'
             for is = 1:numFramesPerModType
                 M = 4;
                 amp = 1 / sqrt(mean(abs(pammod(0:3, 4)).^2));
                 src = randi([0 3],spf/sps,1);
                 if isempty(filterCoeffs)
                  filterCoeffs = rcosdesign(0.35, 4, sps);
                 end
                 if isempty(amp)
                  amp = 1 / sqrt(mean(abs(pammod(0:3, 4)).^2));
                 end
                 % Modulate
                 syms = amp * pammod(src,4);
                  % Pass through independent channels
                 y = channel(syms);
                 
                  % Pulse shape
                 wav = filter(filterCoeffs, 1, upsample(y,sps));
                 
                  % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, spf, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modulationTypes(modType));
                 % Save signal
                 data{idx} = wav;
                 truth(idx) = modulationTypes(modType);
                 idx = idx + 1;
             end
     end
end
end   

  %% Subroutines
function val = randOverInterval(interval)
% Expect interval to be <1x2> with format [minVal maxVal]
val = (interval(2) - interval(1)).*rand + interval(1);
end
    
 
 
    
    
    
% %   dataSrc = getSource(modulationTypes(modType), sps, 2*spf, fs);
%   modulator = getModulator(modulationTypes(modType), sps, fs);
%   if contains(char(modulationTypes(modType)), {'B-FM','DSB-AM','SSB-AM'})
%     % Analog modulation types use 100 MHz center frequency
%     channel.CenterFrequency = 100e6;
%   else
%     % Digital modulation types use 900 MHz center frequency
%     channel.CenterFrequency = 900e6;
%   end
%   
%   for p=1:numFramesPerModType
%     % Generate random data
%     x = dataSrc();
%     
%     % Modulate
%     y = modulator(x);
%     
%     % Pass through independent channels
%     rxSamples = channel(y);
%    
%     % Remove transients from the beginning, trim to size, and normalize
%     data{idx} = rxSamples;
%     truth(idx) = modulationTypes(modType);
%     idx = idx + 1;
%     
%   end
% end
% end
% 
% %% Helper functions
% 
% function modulator = getModulator(modType, sps, fs)
% %getModulator Modulation function selector
% %   MOD = getModulator(TYPE,SPS,FS) returns a modulator function handle,
% %   MOD, based on the TYPE. SPS is samples per symbol and FS is sampling
% %   rate.
% 
% switch modType
%     case "BPSK"
%         modulator = @(x)bpskModulator(x,sps);
%     case "QPSK"
%         modulator = @(x)qpskModulator(x,sps);
%     case "8PSK"
%         modulator = @(x)psk8Modulator(x,sps);
%     case "2FSK"
%        modulator = @(x)fsk2Modulator(x,sps);
%     case "8FSK"
%        modulator = @(x)fsk8Modulator(x,sps);
%     case "16QAM"
%        modulator = @(x)qam16Modulator(x, fs);
%     case "64QAM"
%         modulator = @(x)qam64Modulator(x, fs);
%     case "GFSK"
%         modulator = @(x)gfskModulator(x, fs);
%     case "CPFSK"
%         modulator = @(x)cpfskModulator(x, fs);
%     case "PAM4"
%         modulator = @(x)pam4Modulator(x, fs);
% %     case "B-FM"
% %         modulator = @(x)bfmModulator(x, fs);
% %     case "DSB-AM"
% %         modulator = @(x)dsbamModulator(x, fs);
% 
% end
% end
% 
% % function src = getSource(modType, sps, spf, fs)
% % %getSource Source selector for modulation types
% % %    SRC=getSource(TYPE,SPS,SPF,FS) returns the data source
% % %    for modulation type, TYPE, with samples per symbol, SPS,
% % %    samples per frame, SPF, and sampling frequency, FS.
% % 
% % switch modType
% %    case {"GFSK","CPFSK","2FSK","BPSK"}
% %     M = 2
% %     src = @()randi([0 M-1],spf/sps,1);
% %   case {"QPSK","PAM4"}
% %     M = 4
% %     src = @()randi([0 M-1],spf/sps,1);
% %   case {"8FSK","8PSK"}
% %     M = 8
% %     src = @()randi([0 M-1],spf/sps,1);
% %   case "16QAM"
% %     M = 16
% %     src = @()randi([0 M-1],spf/sps,1);
% %   case "64QAM"
% %     M = 64
% %     src = @()randi([0 M-1],spf/sps,1);
% % end
% % end
% 
% 
% 
% %% Modulators
% 
% function y = bpskModulator(x,sps)
% %bpskModulator BPSK modulator with pulse shaping
% %   Y=bpskModulator(X,SPS) BPSK modulates input, X, and returns a root-raised
% %   cosine pulse shaped signal, Y. X must be a column vector of values in
% %   the set of [0 1]. Root-raised cosine filter has a roll-off factor of
% %   0.35 and spans four symbols. The output signal, Y, has unit power.
% 
% persistent filterCoeffs
% if isempty(filterCoeffs)
%   filterCoeffs = rcosdesign(0.35, 4, sps);
% end
% % Modulate
% syms = pskmod(x,2);
% % Pulse shape
% y = filter(filterCoeffs, 1, upsample(syms,sps));
% end
% 
% function y = qpskModulator(x,sps)
% %qpskModulator QPSK modulator with pulse shaping
% %   Y=qpskModulator(X,SPS) QPSK modulates input, X, and returns a root-raised
% %   cosine pulse shaped signal, Y. X must be a column vector of values in
% %   the set of [0 3]. Root-raised cosine filter has a roll-off factor of
% %   0.35 and spans four symbols. The output signal, Y, has unit power.
% 
% persistent filterCoeffs
% if isempty(filterCoeffs)
%   filterCoeffs = rcosdesign(0.35, 4, sps);
% end
% % Modulate
% syms = pskmod(x,4,pi/4);
% % Pulse shape
% y = filter(filterCoeffs, 1, upsample(syms,sps));
% end
% 
% function y = pam4Modulator(x,sps)
% %pam4Modulator PAM4 modulator with pulse shaping
% %   Y=qam16Modulator(X,SPS) PAM4 modulates input, X, and returns a root-raised
% %   cosine pulse shaped signal, Y. X must be a column vector of values in
% %   the set of [0 3]. Root-raised cosine filter has a roll-off factor of
% %   0.35 and spans four symbols. The output signal, Y, has unit power.
% 
% persistent filterCoeffs amp
% if isempty(filterCoeffs)
%   filterCoeffs = rcosdesign(0.35, 4, sps);
%   amp = 1 / sqrt(mean(abs(pammod(0:3, 4)).^2));
% end
% % Modulate
% syms = amp * pammod(x,4);
% % Pulse shape
% y = filter(filterCoeffs, 1, upsample(syms,sps));
% end
% 
% function y = gfskModulator(x,sps)
% %gfskModulator GFSK modulator
% %   Y=gfskModulator(X,SPS) GFSK modulates input, X, and returns signal, Y. 
% %   X must be a column vector of values in the set of [0 1]. BT product is
% %   0.35 and modulation index is 1. The output signal, Y, has unit power.
% 
% persistent mod meanM
% if isempty(mod)
%   M = 2;
%   mod = comm.CPMModulator(...
%     'ModulationOrder', M, ...
%     'FrequencyPulse', 'Gaussian', ...
%     'BandwidthTimeProduct', 0.35, ...
%     'ModulationIndex', 1, ...
%     'SamplesPerSymbol', sps);
%   meanM = mean(0:M-1);
% end
% % Modulate
% y = mod(2*(x-meanM));
% 
% 
% end
% 
% function y = cpfskModulator(x,sps)
% %cpfskModulator CPFSK modulator
% %   Y=cpfskModulator(X,SPS) CPFSK modulates input, X, and returns signal, Y. 
% %   X must be a column vector of values in the set of [0 1]. Modulation 
% %   index is 0.5. The output signal, Y, has unit power.
% 
% persistent mod meanM
% if isempty(mod)
%   M = 2;
%   mod = comm.CPFSKModulator(...
%     'ModulationOrder', M, ...
%     'ModulationIndex', 0.5, ...
%     'SamplesPerSymbol', sps);
%   meanM = mean(0:M-1);
% end
% % Modulate
% y = mod(2*(x-meanM));
% 
% 
% end
% 
% function y = fsk2Modulator(x,sps)
% %fsk2Modulator FSK modulator with pulse shaping
% %   Y = fsk2Modulator(X,SPS) FSK modulates the input X, and returns the 
% %   root-raised cosine pulse shaped signal Y. X must be a column vector 
% %   of values in the set [0 1]. The root-raised cosine filter has a 
% %   roll-off factor of 0.35 and spans four symbols. The output signal 
% %   Y has unit power.
% 
% persistent mod 
% if isempty(mod)
%   M = 2;
%   freqsep = 4;  % Frequency separation (Hz)
%   symbolRate = 45;
%   fsamp = sps*symbolRate;
%   
%  mod = comm.FSKModulator(...
%     'ModulationOrder', M, ...
%     'SymbolMapping', 'Binary', ...
%     'FrequencySeparation', freqsep, ...
%     'SymbolRate', 45, ...
%     'SamplesPerSymbol', sps);
% 
% 
% end
% % Modulate
% y = mod(x);
% 
% 
% end
% 
% 
% function y = fsk8Modulator(x,sps)
% %fsk8Modulator FSK modulator with pulse shaping
% %   Y = fsk8Modulator(X,SPS) FSK modulates the input X, and returns the 
% %   root-raised cosine pulse shaped signal Y. X must be a column vector 
% %   of values in the set [0 1]. The root-raised cosine filter has a 
% %   roll-off factor of 0.35 and spans four symbols. The output signal 
% %   Y has unit power.
% 
% persistent mod 
% if isempty(mod)
%   M = 8;  
%   freqsep = 4;  % Frequency separation (Hz)
%   symbolRate = 45;
%   
%   mod = comm.FSKModulator(...
%     'ModulationOrder', M, ...
%     'FrequencySeparation', freqsep, ...
%     'SymbolRate', symbolRate, ...
%     'SamplesPerSymbol', sps);
% 
% release(mod);
% end 
% % Modulate
% y = mod(x);
% 
% 
% end
% 
% function y = psk8Modulator(x,sps)
% %psk8Modulator 8-PSK modulator with pulse shaping
% %   Y = psk8Modulator(X,SPS) 8-PSK modulates the input X, and returns the 
% %   root-raised cosine pulse shaped signal Y. X must be a column vector 
% %   of values in the set [0 7]. The root-raised cosine filter has a 
% %   roll-off factor of 0.35 and spans four symbols. The output signal 
% %   Y has unit power.
% 
% persistent filterCoeffs
% if isempty(filterCoeffs)
%   filterCoeffs = rcosdesign(0.35, 4, sps);
% end
% % Modulate
% syms = pskmod(x,8);
% % Pulse shape
% y = filter(filterCoeffs, 1, upsample(syms,sps));
% end
% 
% 
% function y = qam16Modulator(x,sps)
% %print("Generating qam16Modulator....");
% %qam16Modulator 16-QAM modulator with pulse shaping
% %   Y = qam16Modulator(X,SPS) 16-QAM modulates the input X, and returns the 
% %   root-raised cosine pulse shaped signal Y. X must be a column vector 
% %   of values in the set [0 15]. The root-raised cosine filter has a 
% %   roll-off factor of 0.35 and spans four symbols. The output signal 
% %   Y has unit power.
% 
% persistent filterCoeffs
% if isempty(filterCoeffs)
%   filterCoeffs = rcosdesign(0.35, 4, sps);
% end
% % Modulate and pulse shape
% syms = qammod(x,16,'UnitAveragePower',true);
% % Pulse shape
% y = filter(filterCoeffs, 1, upsample(syms,sps));
% end
% 
% function y = qam64Modulator(x,sps)
% %qam64Modulator 64-QAM modulator with pulse shaping
% %   Y = qam64Modulator(X,SPS) 64-QAM modulates the input X, and returns the 
% %   root-raised cosine pulse shaped signal Y. X must be a column vector 
% %   of values in the set [0 63]. The root-raised cosine filter has a 
% %   roll-off factor of 0.35 and spans four symbols. The output signal 
% %   Y has unit power.
% print("Generating qam64Modulator....");
% persistent filterCoeffs
% if isempty(filterCoeffs)
%   filterCoeffs = rcosdesign(0.35, 4, sps);
% end
% % Modulate
% syms = qammod(x,64,'UnitAveragePower',true);
% % Pulse shape
% y = filter(filterCoeffs, 1, upsample(syms,sps));
% end
% 
% function y = bfmModulator(x,fs)
% %bfmModulator Broadcast FM modulator
% %   Y=bfmModulator(X,FS) broadcast FM modulates input, X, and returns 
% %   signal, Y, at a sample rate of FS. X must be a column vector of 
% %   audio samples at sample rate of FS. Frequency deviation is 75 kHz 
% %   and pre-emphasis filter time constant is 75 microseconds. 
% 
% persistent mod
% if isempty(mod)
%   mod = comm.FMBroadcastModulator(...
%     'AudioSampleRate', fs, ...
%     'SampleRate', fs);
% end
% y = mod(x);
% end
% 
% function y = dsbamModulator(x,fs)
% %dsbamModulator Double sideband AM modulator
% %   Y=dsbamModulator(X,FS) double sideband AM modulates input, X, and 
% %   returns signal, Y, at a sample rate of FS. X must be a column vector of 
% %   audio samples at sample rate of FS. IF frequency is 50 kHz. 
% 
% y = ammod(x,50e3,fs);
% end
% 
% function y = ssbamModulator(x,fs)
% %ssbamModulator Single sideband AM modulator
% %   Y=ssbamModulator(X,FS) single sideband AM modulates input, X, and 
% %   returns signal, Y, at a sample rate of FS. X must be a column vector of 
% %   audio samples at sample rate of FS. IF frequency is 50 kHz. 
% 
% y = ssbmod(x,50e3,fs);
% end

