function helperModClassification_CNN(frameStore)

%% Setup Simulation Parameters
% User defined parameters
Fs = 100e3; % Sampling frequency (Hz)
numFramesPerModType = 20; % Number of signals per modulation type
Ts = 1/Fs; % Sampling period (sec)
sps = 8;
spf = 1024;% changed 512 to 1024
windowLen = 1024;% changed 512 to 1024
 persistent filterCoeffs amp mod meanM;
SNR = 20;
%% Initialization
modTypes = categorical(["LFM","Rect","Barker","Frank","P1","P2","P3","P4","Zadoff-Chu",...
    "BPSK","QPSK","8PSK","PAM4","16QAM","64QAM",...
    "GFSK","CPFSK","2FSK","8FSK"]);
% 
% 
 
rng(12345)
tic
hFreqOffset = comm.PhaseFrequencyOffset(...
    'SampleRate',Fs);
%CombinedChannel
channel = helperModClassTestChannel(...
  'SampleRate', Fs, ...
  'SNR', SNR, ...
  'PathDelays', [0 1.8 3.4] / Fs, ...
  'AveragePathGains', [0 -2 -10], ...
  'KFactor', 4, ...
  'MaximumDopplerShift', 4, ...
  'MaximumClockOffset', 5, ...
  'CenterFrequency', 900e6);% Digital modulation types use 900 MHz center frequency
   transDelay = 30;
   
%% Generate training data

for iM = 1:length(modTypes)
    fprintf('%s - Generating %s frames\n', ...
    datestr(toc/86400,'HH:MM:SS'), modTypes(iM))
    modType = modTypes(iM);
    rangeFc = [Fs/6, Fs/5]; % Center frequency (Hz) range
    rangeN = [512, 1920]; % Number of collected signal samples range
    rangeB = [Fs/20, Fs/16]; % Bandwidth (Hz) range
    sweepDirections = {'Up','Down'};
    
    switch modType
        case 'Rect'
            % Create signal
            hRect = phased.RectangularWaveform(...
                'SampleRate',Fs,...
                'OutputFormat','Samples');
            
            for iS = 1:numFramesPerModType
                %Get randomized parameters
                Ncc = round(randOverInterval(rangeN));
                
                % Create waveform
                hRect.PulseWidth = Ncc*Ts;
                hRect.PRF = 1/(Ncc*Ts);
                hRect.NumSamples = 256; % Changed from 1024 to 256
                src = hRect();
                
                
                if isempty(filterCoeffs)
                   filterCoeffs = rcosdesign(0.35, 4, sps);
                end
                % Pulse shaping (Added by Rachana)
                src = filter(filterCoeffs, 1,upsample(src,sps));
                release(channel);
                wav = channel(src);
                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));
                
                release(hRect);
                release(hFreqOffset);
            end
            clear hRect;
        case 'LFM'
            % Create signal
            hLfm = phased.LinearFMWaveform(...
                'SampleRate',Fs,...
                'OutputFormat','Samples');
            
            for iS = 1:numFramesPerModType
                %Get randomized parameters
                B = randOverInterval(rangeB);
                Ncc = round(randOverInterval(rangeN));
                
                % Generate LFM
                hLfm.SweepBandwidth = B;
                hLfm.PulseWidth = Ncc*Ts;
                hLfm.NumSamples = 256;
                hLfm.PRF = 1/(Ncc*Ts);
                hLfm.SweepDirection = sweepDirections{randi(2)};
                src = hLfm();
                
                if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                end
                % Pulse shaping (Added by Rachana)
                wav = filter(filterCoeffs, 1, upsample(src,sps));
                release(channel);
                tx = channel(wav);
                
                
                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(tx, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));
                 
                
                release(hLfm);
                release(hFreqOffset);
            end
            clear hLfm;
        case 'Barker'
            rangeNChip = [3,4,5,7,11]; % Number of chips
            rangeNcc = [1,5]; % Cycles per phase code
            
            % Create signal and update SNR
            hPhase = phased.PhaseCodedWaveform(...
                'SampleRate',Fs,...
                'Code',string(modType),...
                'OutputFormat','Samples');
            
            for iS = 1:numFramesPerModType
                %Get randomized parameters
                Fc = randOverInterval(rangeFc) ;
                N = rangeNChip(randi(length(rangeNChip),1));
                Ncc = rangeNcc(randi(length(rangeNcc),1));

                
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*Fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhase.ChipWidth = chipWidth;
                hPhase.NumChips = N;
                hPhase.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhase.NumSamples = 256;
                wav = hPhase();
                
                if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                end
                % Pulse shaping (Added by Rachana)
                wav = filter(filterCoeffs, 1, upsample(wav,sps));

                wav = channel(wav);
                
                
                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));
                
                release(hPhase);
                release(hFreqOffset);
            end
            clear hPhase;
            
            case 'Frank'
            rangeNChip = 4; % Number of chips
            rangeNcc = [1,5]; % Cycles per phase code
            
            % Create signal and update SNR
            hPhaseFrank = phased.PhaseCodedWaveform(...
                'SampleRate',Fs,...
                'Code',string(modType),...
                'OutputFormat','Samples');
            
            for iS = 1:numFramesPerModType
                %Get randomized parameters
                Fc = randOverInterval(rangeFc) ;
                
                N = rangeNChip(randi(length(rangeNChip),1));
                Ncc = rangeNcc(randi(length(rangeNcc),1));
                
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*Fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhaseFrank.ChipWidth = chipWidth;
                hPhaseFrank.NumChips = N;
                hPhaseFrank.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhaseFrank.NumSamples = 256;
                src = hPhaseFrank();
                
                if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                end
                % Pulse shaping (Added by Rachana)
                wav = filter(filterCoeffs, 1, upsample(src,sps));
                tx = channel(wav);
                
                
                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(tx, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));

                release(hPhaseFrank);
                release(hFreqOffset);
            end
            clear hPhase;
            
         case 'P1'
            rangeNChip = 4; % Number of chips
            rangeNcc = [1,5]; % Cycles per phase code
            
            % Create signal and update SNR
            hPhaseP1 = phased.PhaseCodedWaveform(...
                'SampleRate',Fs,...
                'Code',string(modType),...
                'OutputFormat','Samples');
            
            for iS = 1:numFramesPerModType
                %Get randomized parameters
                Fc = randOverInterval(rangeFc) ;
                N = rangeNChip(randi(length(rangeNChip),1));
                Ncc = rangeNcc(randi(length(rangeNcc),1));
                
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*Fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhaseP1.ChipWidth = chipWidth;
                hPhaseP1.NumChips = N;
                hPhaseP1.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhaseP1.NumSamples = 256;
                src = hPhaseP1();
                if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                end
                % Pulse shaping (Added by Rachana)
                wav = filter(filterCoeffs, 1, upsample(src,sps));

                tx = channel(wav);
                
                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(tx, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));

                release(hPhaseP1);
                release(hFreqOffset);
            end
            clear hPhaseP1;
            
            case 'P2'
            rangeNChip = 4; % Number of chips
            rangeNcc = [1,5]; % Cycles per phase code
            
            % Create signal and update SNR
            hPhaseP2 = phased.PhaseCodedWaveform(...
                'SampleRate',Fs,...
                'Code',string(modType),...
                'OutputFormat','Samples');
            
            for iS = 1:numFramesPerModType
                %Get randomized parameters
                Fc = randOverInterval(rangeFc) ;
                N = rangeNChip(randi(length(rangeNChip),1));
                Ncc = rangeNcc(randi(length(rangeNcc),1));
                
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*Fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhaseP2.ChipWidth = chipWidth;
                hPhaseP2.NumChips = N;
                hPhaseP2.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhaseP2.NumSamples = 256;
                src = hPhaseP2();
                if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                end
               % Pulse shaping (Added by Rachana)
                wav = filter(filterCoeffs, 1, upsample(src,sps));
                
                tx = channel(wav);
                
                
                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(tx, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));
                
                release(hPhaseP2);
                release(hFreqOffset);
            end
            clear hPhaseP2;
         case 'P3'
            rangeNChip = 4; % Number of chips
            rangeNcc = [1,5]; % Cycles per phase code
            
            % Create signal and update SNR
            hPhaseP3 = phased.PhaseCodedWaveform(...
                'SampleRate',Fs,...
                'Code',string(modType),...
                'OutputFormat','Samples');
            
            for iS = 1:numFramesPerModType
                %Get randomized parameters
                Fc = randOverInterval(rangeFc) ;
                N = rangeNChip(randi(length(rangeNChip),1));
                Ncc = rangeNcc(randi(length(rangeNcc),1));
                
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*Fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhaseP3.ChipWidth = chipWidth;
                hPhaseP3.NumChips = N;
                hPhaseP3.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhaseP3.NumSamples = 256;
                src = hPhaseP3();
                if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                end
                % Pulse shaping (Added by Rachana)
                wav = filter(filterCoeffs, 1, upsample(src,sps));

                tx = channel(wav);
                
                
                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(tx, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));
                
                release(hPhaseP3);
                release(hFreqOffset);
            end
            clear hPhaseP3;
            
         case 'P4'
            rangeNChip = 4; % Number of chips
            rangeNcc = [1,5]; % Cycles per phase code
            
            % Create signal and update SNR
            hPhaseP4 = phased.PhaseCodedWaveform(...
                'SampleRate',Fs,...
                'Code',string(modType),...
                'OutputFormat','Samples');
            
            for iS = 1:numFramesPerModType
                %Get randomized parameters
                Fc = randOverInterval(rangeFc) ;
                N = rangeNChip(randi(length(rangeNChip),1));
                Ncc = rangeNcc(randi(length(rangeNcc),1));
                
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*Fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhaseP4.ChipWidth = chipWidth;
                hPhaseP4.NumChips = N;
                hPhaseP4.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhaseP4.NumSamples = 256;
                src = hPhaseP4();
                if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                end
                % Pulse shaping (Added by Rachana)
                wav = filter(filterCoeffs, 1, upsample(src,sps));

                tx = channel(wav);
                
                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(tx, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));
                
                release(hPhaseP4);
                release(hFreqOffset);
            end
            clear hPhaseP4;
        
        case 'Zadoff-Chu'
            rangeNChip = 4; % Number of chips
            rangeNcc = [1,5]; % Cycles per phase code
            
            % Create signal and update SNR
            hPhaseZadoffChu = phased.PhaseCodedWaveform(...
                'SampleRate',Fs,...
                'Code',string(modType),...
                'OutputFormat','Samples');
            
            for iS = 1:numFramesPerModType
                %Get randomized parameters
                Fc = randOverInterval(rangeFc) ;
                N = rangeNChip(randi(length(rangeNChip),1));
                Ncc = rangeNcc(randi(length(rangeNcc),1));
                
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*Fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhaseZadoffChu.ChipWidth = chipWidth;
                hPhaseZadoffChu.NumChips = N;
                hPhaseZadoffChu.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhaseZadoffChu.NumSamples = 256;
                src = hPhaseZadoffChu();
                if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                end
                % Pulse shaping (Added by Rachana)
                wav = filter(filterCoeffs, 1, upsample(src,sps));

                tx = channel(wav);
                
                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(tx, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));

                release(hPhaseZadoffChu);
                release(hFreqOffset);
            end
            clear hPhaseZadoffChu;
          
         case 'BPSK'
             for is = 1:numFramesPerModType
                 M = 2;
                 src = randi([0 M-1],spf*2/sps,1);
                 % bpskModulator BPSK modulator with pulse shaping
                 % Y=bpskModulator(X,SPS) BPSK modulates input, X, and returns a root-raised
                 % cosine pulse shaped signal, Y. X must be a column vector of values in
                 % the set of [0 1]. Root-raised cosine filter has a roll-off factor of
                 % 0.35 and spans four symbols. The output signal, Y, has unit power.
                
                 % Modulate
                 syms = pskmod(src,2);
                  if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                  end
                 
                 % Pulse shape
                 y = filter(filterCoeffs, 1, upsample(syms,sps));
                release(channel);
                % Add multipath offset
                wav = channel(y);
                 
                 % Remove transients from the beginning, trim to size, and normalize
               frame = helperModClassFrameGenerator(wav, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));
                 
                release(hFreqOffset);
             end
         case 'QPSK'
             for is = 1:numFramesPerModType
                 M = 4;
                 src = randi([0 M-1],spf*2/sps,1);
                 
                 % Modulate
                 syms = pskmod(src,4,pi/4);
                 if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                  end
                 
                 % Pulse shape
                 y = filter(filterCoeffs, 1, upsample(syms,sps));
                 
                 % Pass through independent channels
                  wav = channel(y);
                 
                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));
                 
                 release(hFreqOffset)
             end
         case '8PSK'
             for is = 1:numFramesPerModType
                 M = 8;
                 src = randi([0 M-1],spf*2/sps,1);
                 if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                 end
                 % Modulate
                 syms = pskmod(src,8);
                 if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                  end
                 
                 % Pulse shape
                 y = filter(filterCoeffs, 1, upsample(syms,sps));
                 % Pass through independent channels
                  wav = channel(y);

                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));

                 release(hFreqOffset);
             end
         case '2FSK'
             for is = 1:numFramesPerModType
                 M = 2;
                 freqsep = 8;  % Frequency separation (Hz)
                 symbolRate = 45;
%                  fsamp = sps*symbolRate;
                 src = randi([0 M-1],spf*2/sps,1);
                 
                 % Modulate
%                   syms = fskmod(src,M,freqsep,sps,fs);
                 fmod = comm.FSKModulator(...
                       'ModulationOrder', M, ...
                       'FrequencySeparation', freqsep, ...
                       'SymbolRate', symbolRate, ...
                       'SamplesPerSymbol', sps);
                 syms = fmod(src);

                 % Pass through independent channels
                 wav = channel(syms);

                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));

                 release(fmod);
                 release(hFreqOffset);
             end
              clear fmod;
         case '8FSK'
             for is = 1:numFramesPerModType
                 M = 8;
                 src = randi([0 M-1],spf*2/sps,1);
                 
                 % Modulate
                 syms = fskmod(src,M,freqsep,sps,Fs);
                  

                 % Pass through independent channels
                 wav = channel(syms);

                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));
                release(hFreqOffset);
             end
         case 'GFSK'
             for is = 1:numFramesPerModType
                 M = 2;
                 src = randi([0 1],spf*2/sps,1);
                 %   Y=gfskModulator(X,SPS) GFSK modulates input, X, and returns signal, Y. 
                 %   X must be a column vector of values in the set of [0 1]. BT product is
                 %   0.35 and modulation index is 1. The output signal, Y, has unit power.

                 if isempty(mod)
                 gmod = comm.CPMModulator(...
                    'ModulationOrder', M, ...
                    'FrequencyPulse', 'Gaussian', ...
                    'BandwidthTimeProduct', 0.35, ...
                    'ModulationIndex', 1, ...
                    'SamplesPerSymbol', sps);
                 meanM = mean(0:M-1);
                 end
                 % Modulate
                 syms = gmod(2*(src-meanM));
                 
                 % Pass through independent channels
                 wav = channel(syms);
                 
                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));

                 release(gmod);
                 release(hFreqOffset);
             end
             clear gmod;
         case '16QAM'
             for is = 1:numFramesPerModType
                 M = 16;
                 src = randi([0 M-1],spf*2/sps,1);
                 % Modulate and pulse shape
                 syms = qammod(src,16,'UnitAveragePower',true);
                  if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                  end
                 
                 % Pulse shape
                 y = filter(filterCoeffs, 1, upsample(syms,sps));
                 % Pass through independent channels
                 wav = channel(y);
                 
                 % Remove transients from the beginning, trim to size, and normalize
                 frame = helperModClassFrameGenerator(wav, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));

                 release(hFreqOffset);
             end
         case '64QAM'
             for is = 1:numFramesPerModType
                 M = 64;
                 src = randi([0 M-1],spf*2/sps,1);
                 % Modulate and pulse shape
                 syms = qammod(src,64,'UnitAveragePower',true);
                 
                 if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                  end
                 
                 % Pulse shape
                 y = filter(filterCoeffs, 1, upsample(syms,sps));
                 % Pass through independent channels
                wav = channel(y);
              
                % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));

                 release(hFreqOffset);
             end
         case 'CPFSK'
             for is = 1:numFramesPerModType
                 M = 2;
%                    Fc = 900e6 ; %Fc = randOverInterval(rangeFc) ;
                 src = randi([0 1],spf*2/sps,1);
                  if isempty(mod)
                   mod = comm.CPFSKModulator(...
                         'ModulationOrder', M, ...
                         'ModulationIndex', 0.5, ...
                         'SamplesPerSymbol', sps);
                   meanM = mean(0:M-1);
                  end
                 % Modulate
                 syms = mod(2*(src-meanM));

%                  Pass through independent channels
                 wav = channel(syms);

                  % Remove transients from the beginning, trim to size, and normalize
                frame = helperModClassFrameGenerator(wav, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));

                 release(hFreqOffset);
                 release(mod);
             end
             clear mod;
         case 'PAM4'
             for is = 1:numFramesPerModType
                
                 amp = 1 / sqrt(mean(abs(pammod(0:3, 4)).^2));
                 src = randi([0 3],spf*2/sps,1);
                 
                 if isempty(amp)
                  amp = 1 / sqrt(mean(abs(pammod(0:3, 4)).^2));
                 end
                 % Modulate
                 syms = amp * pammod(src,4);
                  if isempty(filterCoeffs)
                    filterCoeffs = rcosdesign(0.35, 4, sps);
                  end
                 
                 % Pulse shape
                 y = filter(filterCoeffs, 1, upsample(syms,sps));
                 % Pass through independent channels
                 wav = channel(y);
                 
                  % Remove transients from the beginning, trim to size, and normalize
                  frame = helperModClassFrameGenerator(wav, windowLen, spf, transDelay, sps);
    
               % Add to frame store
                add(frameStore, frame, modTypes(modType));
                 release(hFreqOffset);
             end
          otherwise
            error('Modulation type not recognized.');
    end
    
end
        
end



%% Subroutines
function val = randOverInterval(interval)
% Expect interval to be <1x2> with format [minVal maxVal]
val = (interval(2) - interval(1)).*rand + interval(1);
end

