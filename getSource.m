function src = getSource(modType, sps, spf, fs)
%getSource Source selector for modulation types
%    SRC = getSource(TYPE,SPS,SPF,fs) returns the data source
%    for the modulation type TYPE, with the number of samples 
%    per symbol SPS, the number of samples per frame SPF, and 
%    the sampling frequency fs.

switch modType
  case {"BPSK","2FSK","GFSK","CPFSK"}
    M = 2;
    src = @()randi([0 M-1],spf/sps,1);
  case {"QPSK","PAM4"}
    M = 4;
    src = @()randi([0 M-1],spf/sps,1);
  case {"8PSK","8FSK"}
    M = 8;
    src = @()randi([0 M-1],spf/sps,1);
  case "16QAM"
    M = 16;
    src = @()randi([0 M-1],spf/sps,1);
  case "64QAM"
    M = 64;
    src = @()randi([0 M-1],spf/sps,1);
  case {"B-FM","DSB-AM","SSB-AM"}
    src = @()getAudio(spf,fs);
  case "LFM"
%     rangeN = [1024, 2432]; % Number of collected signal samples range
    rangeN = [512, 1920];
    rangeB = [fs/20, fs/16]; % Bandwidth (Hz) range
    sweepDirections = {'Up','Down'};
    Ts = 1/fs;
    
    %Get randomized parameters
    B = randOverInterval(rangeB);
    Ncc = round(randOverInterval(rangeN));
    hLfm = phased.LinearFMWaveform('SampleRate',fs,'OutputFormat','Samples');           
    % Generate LFM
    hLfm.SweepBandwidth = B;
    hLfm.PulseWidth = Ncc*Ts;
    hLfm.NumSamples = 256;
    hLfm.PRF = 1/(Ncc*Ts);
    hLfm.SweepDirection = sweepDirections{randi(2)};
    src = hLfm();  
   
%    filter = phased.MatchedFilter( ...
%     'Coefficients',getMatchedFilter(hLfm),...
%     'SampleRate', fs,...
%     'SpectrumWindow','None');
%    src = filter(wav);
   
   case 'Rect'
    % Create signal
    hRect = phased.RectangularWaveform(...
             'SampleRate',fs,...
             'OutputFormat','Samples');
            
    %Get randomized parameters
    rangeN = [512, 1920]; % Number of collected signal samples range
    Ts = 1/fs;
    Ncc = round(randOverInterval(rangeN));
                
   % Create waveform
    hRect.PulseWidth = Ncc*Ts;
    hRect.PRF = 1/(Ncc*Ts);
    hRect.NumSamples = 256;
    src = hRect();
                
%     filter = phased.MatchedFilter( ...
%                'Coefficients',getMatchedFilter(hRect),...
%                'SampleRate', fs,...
%                'SpectrumWindow','None');
%                 src = filter(wav);
                
   case 'Barker'
            rangeNChip = [3,4,5,7,11]; % Number of chips
            rangeNcc = [1,5]; % Cycles per phase code
            rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
            Ts = 1/fs;
            % Create signal and update SNR
            hPhaseBarker = phased.PhaseCodedWaveform(...
                'SampleRate',fs,...
                'Code',string(modType),...
                'OutputFormat','Samples');
            
%                 Get randomized parameters
                Fc = randOverInterval(rangeFc);
                N = rangeNChip(randi(length(rangeNChip),1));
                Ncc = rangeNcc(randi(length(rangeNcc),1));
                
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhaseBarker.ChipWidth = chipWidth;
                hPhaseBarker.NumChips = N;
                hPhaseBarker.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhaseBarker.NumSamples = 256;
                src = hPhaseBarker();
                
%                 filter = phased.MatchedFilter( ...
%                           'Coefficients',getMatchedFilter(hPhaseBarker),...
%                           'SampleRate', fs,...
%                           'SpectrumWindow','None');
%                 src = filter(wav);
                
      case 'Frank'
            rangeNChip = 4; % Number of chips
            rangeNcc = [1,5]; % Cycles per phase code
            rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
            Ts = 1/fs;
            % Create signal and update SNR
            hPhaseFrank = phased.PhaseCodedWaveform(...
                'SampleRate',fs,...
                'Code',string(modType),...
                'OutputFormat','Samples');
            
            
%               Get randomized parameters
                
                Fc = randOverInterval(rangeFc);
                N = rangeNChip(randi(length(rangeNChip),1));
                Ncc = rangeNcc(randi(length(rangeNcc),1));
                
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhaseFrank.ChipWidth = chipWidth;
                hPhaseFrank.NumChips = N;
                hPhaseFrank.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhaseFrank.NumSamples = 256;
                src = hPhaseFrank();
                
%                 filter = phased.MatchedFilter( ...
%                           'Coefficients',getMatchedFilter(hPhaseFrank),...
%                           'SampleRate', fs,...
%                           'SpectrumWindow','None');
%                 src = filter(wav);
            
      case 'P1'
            rangeNChip = 4; % Number of chips
            rangeNcc = [1,5]; % Cycles per phase code
            
            % Create signal and update SNR
            hPhaseP1 = phased.PhaseCodedWaveform(...
                'SampleRate',fs,...
                'Code',string(modType),...
                'OutputFormat','Samples');
           rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
            Ts = 1/fs; 
           
          
                %Get randomized parameters
                Fc = randOverInterval(rangeFc);
                N = rangeNChip(randi(length(rangeNChip),1));
                Ncc = rangeNcc(randi(length(rangeNcc),1));

                          
           % Create signal and update SNR
            chipWidth = Ncc/Fc;
            chipWidthSamples = round(chipWidth*fs)-1; % This must be an integer!
            chipWidth = chipWidthSamples*Ts;
            hPhaseP1.ChipWidth = chipWidth;
            hPhaseP1.NumChips = N;
            hPhaseP1.PRF = 1/((chipWidthSamples*N+1)*Ts);
            hPhaseP1.NumSamples = 256;
            src = hPhaseP1();
            
%             filter = phased.MatchedFilter( ...
%                           'Coefficients',getMatchedFilter(hPhaseP1),...
%                           'SampleRate', fs,...
%                           'SpectrumWindow','None');
%             src = filter(wav);
            
            
      case 'P2'
            rangeNChip = 4; % Number of chips
            rangeNcc = [1,5]; % Cycles per phase code
            rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
            Ts = 1/fs; 
            % Create signal and update SNR
            hPhaseP2 = phased.PhaseCodedWaveform(...
                'SampleRate',fs,...
                'Code',string(modType),...
                'OutputFormat','Samples');
            
           
                %Get randomized parameters
               
                Fc = randOverInterval(rangeFc);
                N = rangeNChip(randi(length(rangeNChip),1));
                Ncc = rangeNcc(randi(length(rangeNcc),1));

                            
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhaseP2.ChipWidth = chipWidth;
                hPhaseP2.NumChips = N;
                hPhaseP2.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhaseP2.NumSamples = 256;
                src = hPhaseP2();
                
%                 filter = phased.MatchedFilter( ...
%                           'Coefficients',getMatchedFilter(hPhaseP2),...
%                           'SampleRate', fs,...
%                           'SpectrumWindow','None');
%                 src = filter(wav);
               
     case 'P3'
            rangeNChip = 4; % Number of chips
            rangeNcc = [1,5]; % Cycles per phase code
            
            % Create signal and update SNR
            hPhaseP3 = phased.PhaseCodedWaveform(...
                'SampleRate',fs,...
                'Code',string(modType),...
                'OutputFormat','Samples');
            rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
            Ts = 1/fs; 
           
               
                %Get randomized parameters
                Fc = randOverInterval(rangeFc);
                N = rangeNChip(randi(length(rangeNChip),1));
                Ncc = rangeNcc(randi(length(rangeNcc),1));

           
                
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhaseP3.ChipWidth = chipWidth;
                hPhaseP3.NumChips = N;
                hPhaseP3.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhaseP3.NumSamples = 256;
                src = hPhaseP3();
                
%                 filter = phased.MatchedFilter( ...
%                           'Coefficients',getMatchedFilter(hPhaseP3),...
%                           'SampleRate', fs,...
%                           'SpectrumWindow','None');
%                 src = filter(wav);
            
            case 'P4'
            rangeNChip = 4; % Number of chips
            rangeNcc = [1,5]; % Cycles per phase code
            rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
            Ts = 1/fs; 
            % Create signal and update SNR
            hPhaseP4 = phased.PhaseCodedWaveform(...
                'SampleRate',fs,...
                'Code',string(modType),...
                'OutputFormat','Samples');
            
           
               
                %Get randomized parameters
                Fc = randOverInterval(rangeFc);
                N = rangeNChip(randi(length(rangeNChip),1));
                Ncc = rangeNcc(randi(length(rangeNcc),1));

         
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhaseP4.ChipWidth = chipWidth;
                hPhaseP4.NumChips = N;
                hPhaseP4.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhaseP4.NumSamples = 256;
                src = hPhaseP4();
                
%                 filter = phased.MatchedFilter( ...
%                           'Coefficients',getMatchedFilter(hPhaseP4),...
%                           'SampleRate', fs,...
%                           'SpectrumWindow','None');
%                src = filter(wav);
            
     case 'Zadoff-Chu'
            rangeNChip = 4; % Number of chips
            rangeNcc = [1,5]; % Cycles per phase code
            rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
            Ts = 1/fs; 
            % Create signal and update SNR
            hPhaseZadoffChu = phased.PhaseCodedWaveform(...
                'SampleRate',fs,...
                'Code',string(modType),...
                'OutputFormat','Samples');
            
            
                
                %Get randomized parameters
                Fc = randOverInterval(rangeFc);
                N = rangeNChip(randi(length(rangeNChip),1));
                Ncc = rangeNcc(randi(length(rangeNcc),1));

                
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhaseZadoffChu.ChipWidth = chipWidth;
                hPhaseZadoffChu.NumChips = N;
                hPhaseZadoffChu.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhaseZadoffChu.NumSamples = 256;
                src = hPhaseZadoffChu();
                
%                 filter = phased.MatchedFilter( ...
%                           'Coefficients',getMatchedFilter(hPhaseZadoffChu),...
%                           'SampleRate', fs,...
%                           'SpectrumWindow','None');
%             src = filter(wav);
        otherwise
            error('Modulation type not recognized.');
end


%% Subroutines
function val = randOverInterval(interval)
% Expect interval to be <1x2> with format [minVal maxVal]
val = (interval(2) - interval(1)).*rand + interval(1);
end
end
