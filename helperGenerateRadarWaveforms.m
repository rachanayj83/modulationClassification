function [data, truth] = helperGenerateRadarWaveforms()

%% Setup Simulation Parameters
% User defined parameters
Fs = 1e8; % Sampling frequency (Hz)
nSignalsPerMod = 3000; % Number of signals per modulation type
Ts = 1/Fs; % Sampling period (sec)

%% Initialization
modTypes = categorical(["LFM","Rect","Barker"]);

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
for iM = 1:length(modTypes)
    modType = modTypes(iM);
    rangeFc = [Fs/6, Fs/5]; % Center frequency (Hz) range
    rangeN = [512, 1920]; % Number of collected signal samples range
    snrVector = -6:5:30; % Range of Signal-to-Noise Ratios (SNRs) for training (dB)
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
                SNR = snrVector(randi(length(snrVector),1));
                
                % Create waveform
                hRect.PulseWidth = Ncc*Ts;
                hRect.PRF = 1/(Ncc*Ts);
                hRect.NumSamples = 1024;
                wav = hRect();
                
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
                SNR = snrVector(randi(length(snrVector),1));
                
                % Generate LFM
                hLfm.SweepBandwidth = B;
                hLfm.PulseWidth = Ncc*Ts;
                hLfm.NumSamples = 1024;
                hLfm.PRF = 1/(Ncc*Ts);
                hLfm.SweepDirection = sweepDirections{randi(2)};
                wav = hLfm();
                
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
                SNR = snrVector(randi(length(snrVector),1));
                
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*Fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhase.ChipWidth = chipWidth;
                hPhase.NumChips = N;
                hPhase.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhase.NumSamples = 1024;
                wav = hPhase();
                
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
                release(hPhase);
                release(hFreqOffset);
            end
            clear hPhase;
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