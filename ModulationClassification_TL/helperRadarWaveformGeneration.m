function [data, truth] = helperRadarWaveformGeneration()

%% Setup Simulation Parameters
% User defined parameters
Fs = 100e6; % Sampling frequency (Hz)
nSignalsPerMod = 2; % Number of signals per modulation type
Ts = 1/Fs; % Sampling period (sec)

%% Initialization
modTypes = categorical(["LFM","Rect","Barker","Frank","P1","P2","P3","P4","Zadoff-Chu"]);

rng(0)

multipathChannel = comm.RicianChannel(...
    'SampleRate', Fs, ...
    'PathDelays', [0 1.8 3.4]/Fs, ...
    'AveragePathGains', [0 -2 -10], ...
    'KFactor', 4, ...
    'MaximumDopplerShift', 4);

hFreqOffset = comm.PhaseFrequencyOffset(...
    'SampleRate',Fs);

%% Generate training data
idxW = 1;
sps = 8;
filterCoeffs = rcosdesign(0.35, 4, sps);
for iM = 1:length(modTypes)
    modType = modTypes(iM);
    rangeFc = [Fs/6, Fs/5]; % Center frequency (Hz) range
    rangeN = [512, 1920]; % Number of collected signal samples range
%     snrvector = -20:5:20; % Range of Signal-to-Noise Ratios (SNRs) for training (dB)
    SNR = 20;
    rangeB = [Fs/20, Fs/16]; % Bandwidth (Hz) range
    sweepDirections = {'Up','Down'};
    
    switch modType
        case 'Rect'
            % Create signal
            hRect = phased.RectangularWaveform(...
                'SampleRate',Fs,...
                'OutputFormat','Samples');
            
            for iS = 1:nSignalsPerMod
                %Get randomized parameters
                Fc = randOverInterval(rangeFc);
                Ncc = round(randOverInterval(rangeN));
                
                % Create waveform
                hRect.PulseWidth = Ncc*Ts;
                hRect.PRF = 1/(Ncc*Ts);
                hRect.NumSamples = 256;
                wav = hRect();
                
                % Pulse shape
                wav = filter(filterCoeffs, 1, upsample(wav,sps));
                
                % Adjust SNR
                wav = awgn(wav,SNR);
                
                % Add frequency offset
                hFreqOffset.FrequencyOffset = Fc;
                wav = hFreqOffset(wav); % Frequency shift
                
                % Add multipath offset
                wav = multipathChannel(wav);
                
                % Save signal
                data{idxW} = wav;
                truth(idxW) = modType;
                
                idxW = idxW + 1;
                release(hRect);
                release(hFreqOffset);
            end
            clear hRect;
        case 'LFM'
            % Create signal
            hLfm = phased.LinearFMWaveform(...
                'SampleRate',Fs,...
                'OutputFormat','Samples');
            
            for iS = 1:nSignalsPerMod
                %Get randomized parameters
                Fc = randOverInterval(rangeFc);
                B = randOverInterval(rangeB);
                Ncc = round(randOverInterval(rangeN));
                
                % Generate LFM
                hLfm.SweepBandwidth = B;
                hLfm.PulseWidth = Ncc*Ts;
                hLfm.NumSamples = 256;
                hLfm.PRF = 1/(Ncc*Ts);
                hLfm.SweepDirection = sweepDirections{randi(2)};
                wav = hLfm();
                
                % Pulse shape
                wav = filter(filterCoeffs, 1, upsample(wav,sps));
                
                % Adjust SNR
                wav = awgn(wav,SNR);
                
                % Add frequency offset
                hFreqOffset.FrequencyOffset = Fc;
                wav = hFreqOffset(wav); % Frequency shift
                
                % Add multipath offset
                wav = multipathChannel(wav);
                
                % Save signal
                data{idxW} = wav;
                truth(idxW) = modType;
                
                idxW = idxW + 1;
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
            
            for iS = 1:nSignalsPerMod
                %Get randomized parameters
                Fc = randOverInterval(rangeFc);
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
                
                % Adjust SNR
                wav = awgn(wav,SNR);
                
                % Pulse shape
                wav = filter(filterCoeffs, 1, upsample(wav,sps));
                
                % Add frequency offset
                hFreqOffset.FrequencyOffset = Fc;
                wav = hFreqOffset(wav); % Frequency shift
                
                % Add multipath offset
                wav = multipathChannel(wav);
                
                % Save signal
                data{idxW} = wav;
                truth(idxW) = modType;
                
                idxW = idxW + 1;
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
            
            for iS = 1:nSignalsPerMod
                %Get randomized parameters
                Fc = randOverInterval(rangeFc);
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
                wav = hPhaseFrank();
                
                % Pulse shape
                wav = filter(filterCoeffs, 1, upsample(wav,sps));
                
                % Adjust SNR
                wav = awgn(wav,SNR);
                
                % Add frequency offset
                hFreqOffset.FrequencyOffset = Fc;
                wav = hFreqOffset(wav); % Frequency shift
                
                % Add multipath offset
                wav = multipathChannel(wav);
                
                % Save signal
                data{idxW} = wav;
                truth(idxW) = modType;
                
                idxW = idxW + 1;
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
            
            for iS = 1:nSignalsPerMod
                %Get randomized parameters
                Fc = randOverInterval(rangeFc);
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
                wav = hPhaseP1();
                
                % Pulse shape
                wav = filter(filterCoeffs, 1, upsample(wav,sps));
                
                % Adjust SNR
                wav = awgn(wav,SNR);
                
                % Add frequency offset
                hFreqOffset.FrequencyOffset = Fc;
                wav = hFreqOffset(wav); % Frequency shift
                
                % Add multipath offset
                wav = multipathChannel(wav);
                
                % Save signal
                data{idxW} = wav;
                truth(idxW) = modType;
                
                idxW = idxW + 1;
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
            
            for iS = 1:nSignalsPerMod
                %Get randomized parameters
                Fc = randOverInterval(rangeFc);
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
                wav = hPhaseP2();
                
                % Pulse shape
                wav = filter(filterCoeffs, 1, upsample(wav,sps));
                
                % Adjust SNR
                wav = awgn(wav,SNR);
                
                % Add frequency offset
                hFreqOffset.FrequencyOffset = Fc;
                wav = hFreqOffset(wav); % Frequency shift
                
                % Add multipath offset
                wav = multipathChannel(wav);
                
                % Save signal
                data{idxW} = wav;
                truth(idxW) = modType;
                
                idxW = idxW + 1;
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
            
            for iS = 1:nSignalsPerMod
                %Get randomized parameters
                Fc = randOverInterval(rangeFc);
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
                wav = hPhaseP3();
                
                % Pulse shape
                wav = filter(filterCoeffs, 1, upsample(wav,sps));
                
                % Adjust SNR
                wav = awgn(wav,SNR);
                
                % Add frequency offset
                hFreqOffset.FrequencyOffset = Fc;
                wav = hFreqOffset(wav); % Frequency shift
                
                % Add multipath offset
                wav = multipathChannel(wav);
                
                % Save signal
                data{idxW} = wav;
                truth(idxW) = modType;
                
                idxW = idxW + 1;
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
            
            for iS = 1:nSignalsPerMod
                %Get randomized parameters
                Fc = randOverInterval(rangeFc);
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
                wav = hPhaseP4();
                
                % Pulse shape
                wav = filter(filterCoeffs, 1, upsample(wav,sps));
                
                % Adjust SNR
                wav = awgn(wav,SNR);
                
                % Add frequency offset
                hFreqOffset.FrequencyOffset = Fc;
                wav = hFreqOffset(wav); % Frequency shift
                
                % Add multipath offset
                wav = multipathChannel(wav);
                
                % Save signal
                data{idxW} = wav;
                truth(idxW) = modType;
                
                idxW = idxW + 1;
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
            
            for iS = 1:nSignalsPerMod
                %Get randomized parameters
                Fc = randOverInterval(rangeFc);
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
                wav = hPhaseZadoffChu();
                
                % Pulse shape
                wav = filter(filterCoeffs, 1, upsample(wav,sps));
                
                % Adjust SNR
                wav = awgn(wav,SNR);
                
                % Add frequency offset
                hFreqOffset.FrequencyOffset = Fc;
                wav = hFreqOffset(wav); % Frequency shift
                
                % Add multipath offset
                wav = multipathChannel(wav);
                
                % Save signal
                data{idxW} = wav;
                truth(idxW) = modType;
                
                idxW = idxW + 1;
                release(hPhaseZadoffChu);
                release(hFreqOffset);
            end
            clear hPhaseZadoffChu;
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

