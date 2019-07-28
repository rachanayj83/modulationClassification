modulationTypes = categorical(["BPSK", "QPSK", "8PSK", ...
  "2FSK","8FSK","GFSK", "CPFSK", ...
  "16QAM", "64QAM", "PAM4","LFM",...
  "Rect","Barker","Frank","P1","P2",...
"P3","P4","Zadoff-Chu"]); 
%   "B-FM", "DSB-AM", "SSB-AM"]);

SNR  = -20;
fs = 1e6;

%%
 load trainedNet20SNR_1a
%% BPSK
% Set the random number generator to a known state to be able to regenerate
% the same frames every time the simulation is run
rng(123456)
% Random bits
d = randi([0 1],1024,1);
% BPSK modulation
syms = pskmod(d,2);
% Square-root raised cosine filter
RCfilterCoeffs = rcosdesign(0.35,4,8);
tx = filter(RCfilterCoeffs,1,upsample(syms,8));

% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)
% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction1,score1] = classify(trainedNet20SNR_1a,unknownFrames);

prediction1

plotScores(score1,modulationTypes)

 %% PAM4

% Random bits
d = randi([0 3], 1024, 1);
% PAM4 modulation
syms = pammod(d,4);
% Square-root raised cosine filter
filterCoeffs = rcosdesign(0.35, 4, 8);
tx = filter(filterCoeffs, 1, upsample(syms,8));
% Channel
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-2 2],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)



% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[estimate2,score2] = classify(trainedNet20SNR_1a,unknownFrames);
estimate2

plotScores(score2,modulationTypes)

%% QPSK
d = randi([0 3],1024,1);
% QPSK modulation
syms = pskmod(d,4,pi/4);
% Square-root raised cosine filter
filterCoeffs = rcosdesign(0.35,4,8);
tx = filter(filterCoeffs,1,upsample(syms,8));
% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)


% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction3,score3] = classify(trainedNet20SNR_1a,unknownFrames);


prediction3

plotScores(score3,modulationTypes)


%% 8PSK
d = randi([0 7],1024,1);
% 8PSK modulation
syms = pskmod(d,8);
% Square-root raised cosine filter
filterCoeffs = rcosdesign(0.35,4,8);
tx = filter(filterCoeffs,1,upsample(syms,8));
% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)


% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction4,score4] = classify(trainedNet20SNR_1a,unknownFrames);


prediction4

plotScores(score4,modulationTypes)

%%
d = randi([0 1],1024,1);
% 2FSK modulation
M = 2;        % Modulation order
freqsep = 4;  % Frequency separation (Hz)
symbolRate = 45;
sps = 8;
%fsamp = sps*symbolRate;
  
 mod = comm.FSKModulator(...
    'ModulationOrder', M, ...
    'SymbolMapping', 'Binary', ...
    'FrequencySeparation', freqsep, ...
    'SymbolRate', symbolRate, ...
    'SamplesPerSymbol', sps);

% Modulate
tx = mod(d);
% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)


% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction5,score5] = classify(trainedNet20SNR_1a,unknownFrames);


prediction5


plotScores(score5,modulationTypes)

%%
d = randi([0 7],1024,1);
% 8FSK modulation
M = 8;        % Modulation order
freqsep = 4;  % Frequency separation (Hz)
symbolRate = 45;
sps = 8;
%fsamp = sps*symbolRate;
  
 mod = comm.FSKModulator(...
    'ModulationOrder', M, ...
    'SymbolMapping', 'Binary', ...
    'FrequencySeparation', freqsep, ...
    'SymbolRate', symbolRate, ...
    'SamplesPerSymbol', sps);

% Modulate
tx = mod(d);
% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)


% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction6,score6] = classify(trainedNet20SNR_1a,unknownFrames);

prediction6

plotScores(score6,modulationTypes)

%%
d = randi([0 1],1024,1);
% GFSK modulation
M = 2;        % Modulation order
freqsep = 4;  % Frequency separation (Hz)
symbolRate = 45;
sps = 8;
%fsamp = sps*symbolRate;
  
  mod = comm.CPMModulator(...
    'ModulationOrder', M, ...
    'FrequencyPulse', 'Gaussian', ...
    'BandwidthTimeProduct', 0.35, ...
    'ModulationIndex', 1, ...
    'SamplesPerSymbol', sps);
  meanM = mean(0:M-1);

% Modulate
tx = mod(2*(d-meanM));

% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)

% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction7,score7] = classify(trainedNet20SNR_1a,unknownFrames);


prediction7

plotScores(score7,modulationTypes)

%%
d = randi([0 1],1024,1);
% CPFSK modulation
 M = 2;
  mod = comm.CPFSKModulator(...
    'ModulationOrder', M, ...
    'ModulationIndex', 0.5, ...
    'SamplesPerSymbol', sps);
  meanM = mean(0:M-1);

% Modulate
tx = mod(2*(d-meanM));

% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)


% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction8,score8] = classify(trainedNet20SNR_1a,unknownFrames);


prediction8

plotScores(score8,modulationTypes)


%% % 16QAM modulation
d = randi([0 15],1024,1);
sps = 8;
filterCoeffs = rcosdesign(0.35, 4, sps);

syms = qammod(d,16,'UnitAveragePower',true);
% Pulse shape
tx = filter(filterCoeffs, 1, upsample(syms,sps));
% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)


% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction9,score9] = classify(trainedNet20SNR_1a,unknownFrames);


prediction9

plotScores(score9,modulationTypes)


%% % 64QAM modulation
d = randi([0 63],1024,1);
sps = 8;
filterCoeffs = rcosdesign(0.35, 4, sps);

syms = qammod(d,64,'UnitAveragePower',true);
% Pulse shape
tx = filter(filterCoeffs, 1, upsample(syms,sps));


% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)


% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction10,score10] = classify(trainedNet20SNR_1a,unknownFrames);

prediction10

plotScores(score10,modulationTypes)

%% LFM

rangeN = [512, 1920]; % Number of collected signal samples range
rangeB = [fs/20, fs/16]; % Bandwidth (Hz) range
sweepDirections = {'Up','Down'};
Ts = 1/fs;
sps = 8;
    
%Get randomized parameters
B = randOverInterval(rangeB);
Ncc = round(randOverInterval(rangeN));
    
hLfm = phased.LinearFMWaveform('SampleRate',fs,'OutputFormat','Samples');           
% Generate LFM
hLfm.SweepBandwidth = B;
hLfm.PulseWidth = Ncc*Ts;
hLfm.NumSamples = 1024;
hLfm.PRF = 1/(Ncc*Ts);
hLfm.SweepDirection = sweepDirections{randi(2)};
src = hLfm(); 

% % Pulse shape
% filter = phased.MatchedFilter( ...
%     'Coefficients',getMatchedFilter(hLfm),...
%     'SampleRate', fs,...
%     'SpectrumWindow','None');
% tx = filter(src); 

% Square-root raised cosine filter
filterCoeffs = rcosdesign(0.35, 4, 8);
tx = filter(filterCoeffs, 1, upsample(src,sps));

   
% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);

rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)

% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction11,score11] = classify(trainedNet20SNR_1a,unknownFrames);


prediction11

plotScores(score11,modulationTypes)

%% Rect

sps = 8;
rangeN = [512, 1920]; % Number of collected signal samples range
rangeB = [fs/20, fs/16]; % Bandwidth (Hz) range
sweepDirections = {'Up','Down'};
Ts = 1/fs;
   
%Get randomized parameters
B = randOverInterval(rangeB);
Ncc = round(randOverInterval(rangeN));

% Create signal
            hRect = phased.RectangularWaveform(...
                'SampleRate',fs,...
                'OutputFormat','Samples');
                
                % Create waveform
                hRect.PulseWidth = Ncc*Ts;
                hRect.PRF = 1/(Ncc*Ts);
                hRect.NumSamples = 512;
                src = hRect();

% % Pulse shape
% filter = phased.MatchedFilter( ...
%     'Coefficients',getMatchedFilter(hRect),...
%     'SpectrumWindow','None');
% tx = filter(src);
filterCoeffs = rcosdesign(0.35, 4, sps);
tx = filter(filterCoeffs, 1, upsample(src,sps));

% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)


% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction12,score12] = classify(trainedNet20SNR_1a,unknownFrames);


prediction12

plotScores(score12,modulationTypes)

%%  'Barker'
sps = 8;
rangeNChip = [3,4,5,7,11]; % Number of chips
rangeNcc = [1,5]; % Cycles per phase code
rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
Ts = 1/fs;
% Create signal and update SNR
hPhaseBarker = phased.PhaseCodedWaveform(...
          'SampleRate',fs,...
          'Code','Barker',...
          'OutputFormat','Samples');
            
%Get randomized parameters
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
hPhaseBarker.NumSamples = 512;
src = hPhaseBarker();
                
% Pulse shape
% filter = phased.MatchedFilter( ...
%     'Coefficients',getMatchedFilter(hPhaseBarker),...
%     'SpectrumWindow','None');
% tx = filter(src); 
filterCoeffs = rcosdesign(0.35, 4, sps);
tx = filter(filterCoeffs, 1, upsample(src,sps));

% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
% scope(tx,rx)


% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction13,score13] = classify(trainedNet20SNR_1a,unknownFrames);

prediction13

% plotScores(score13,modulationTypes)
                
                
%%  'Frank'
rangeNChip = 4; % Number of chips
rangeNcc = [1,5]; % Cycles per phase code
rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
Ts = 1/fs;
% Create signal and update SNR
hPhaseFrank = phased.PhaseCodedWaveform(...
                'SampleRate',fs,...
                'Code','Frank',...
                'OutputFormat','Samples');
            
            
%Get randomized parameters
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
hPhaseFrank.NumSamples = 512;
src = hPhaseFrank();
                
% Pulse shape
% filter = phased.MatchedFilter( ...
%             'Coefficients',getMatchedFilter(hPhaseFrank),...
%             'SampleRate', fs,...
%             'SpectrumWindow','None');
% tx = filter(src);
filterCoeffs = rcosdesign(0.35, 4, sps);
tx = filter(filterCoeffs, 1, upsample(src,sps));
                
% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
% scope(tx,rx)

% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction14,score14] = classify(trainedNet20SNR_1a,unknownFrames);

prediction14

% plotScores(score14,modulationTypes)
            
%%  'P1'
rangeNChip = 4; % Number of chips
rangeNcc = [1,5]; % Cycles per phase code
sps = 8;          
% Create signal and update SNR
hPhaseP1 = phased.PhaseCodedWaveform(...
                'SampleRate',fs,...
                'Code','P1',...
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
hPhaseP1.NumSamples = 1024;
src = hPhaseP1();
% bw = bandwidth(hPhaseP1) % Bandwidth of signal

% Pulse shape
% filter = phased.MatchedFilter( ...
%     'Coefficients',getMatchedFilter(hPhaseP1),...
%     'SpectrumWindow','None');
% tx = filter(src);    
filterCoeffs = rcosdesign(0.35, 4, sps);
tx = filter(filterCoeffs, 1, upsample(src,sps));

% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
% scope(tx,rx)


% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction15,score15] = classify(trainedNet20SNR_1a,unknownFrames);


prediction15

% plotScores(score15,modulationTypes)
            
%% 'P2'
sps = 8;
rangeNChip = 4; % Number of chips
rangeNcc = [1,5]; % Cycles per phase code
rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
Ts = 1/fs; 
 % Create signal and update SNR
hPhaseP2 = phased.PhaseCodedWaveform(...
                'SampleRate',fs,...
                'Code','P2',...
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
hPhaseP2.NumSamples = 1024;
src = hPhaseP2();
% bw = bandwidth(hPhaseP2) % Bandwidth of signal

% Pulse shape
% filter = phased.MatchedFilter( ...
%     'Coefficients',getMatchedFilter(hPhaseP2),...
%     'SampleRate', fs,...
%     'SpectrumWindow','None');
% tx = filter(src);
filterCoeffs = rcosdesign(0.35, 4, sps);
tx = filter(filterCoeffs, 1, upsample(src,sps));

% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(src,rx)


% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction16,score16] = classify(trainedNet20SNR_1a,unknownFrames);

prediction16

plotScores(score16,modulationTypes)
               
%% 'P3'
sps = 8;
rangeNChip = 4; % Number of chips
rangeNcc = [1,5]; % Cycles per phase code
rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
Ts = 1/fs; 
 % Create signal and update SNR
hPhaseP3 = phased.PhaseCodedWaveform(...
                'SampleRate',fs,...
                'Code','P3',...
                'OutputFormat','Samples');
            
           
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
hPhaseP3.NumSamples = 512;
src = hPhaseP3();

% Pulse shape
% filter = phased.MatchedFilter( ...
%     'Coefficients',getMatchedFilter(hPhaseP3),...
%     'SampleRate', fs,...
%     'SpectrumWindow','None');
% tx = filter(src);
filterCoeffs = rcosdesign(0.35, 4, sps);
tx = filter(filterCoeffs, 1, upsample(src,sps));

% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)

% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');

% Classification
[prediction17,score17] = classify(trainedNet20SNR_1a,unknownFrames);

prediction17

plotScores(score17,modulationTypes)
            
%% 'P4'
sps = 8;
rangeNChip = 5; % Number of chips
rangeNcc = [1,5]; % Cycles per phase code
rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
Ts = 1/fs; 
 % Create signal and update SNR
hPhaseP4 = phased.PhaseCodedWaveform(...
                'SampleRate',fs,...
                'Code','P4',...
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
hPhaseP4.NumSamples = 1024;
src = hPhaseP4();

% Pulse shape
% filter = phased.MatchedFilter( ...
%     'Coefficients',getMatchedFilter(hPhaseP4),...
%     'SampleRate', fs,...
%     'SpectrumWindow','None');
% tx = filter(src);  
filterCoeffs = rcosdesign(0.35, 4, sps);
tx = filter(filterCoeffs, 1, upsample(src,sps));

% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);

% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)


% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction18,score18] = classify(trainedNet20SNR_1a,unknownFrames);

prediction18

 plotScores(score18,modulationTypes)
            
 %% 'Zadoff-Chu'
sps = 8; 
fs = 1e6;
% rangeNChip = 4; % Number of chips
% rangeNcc = [1,5]; % Cycles per phase code
% rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
Ts = 1/fs; 
% Create signal and update SNR
hPhaseZadoffChu = phased.PhaseCodedWaveform(...
                'SampleRate',fs,...
                'Code','Zadoff-Chu',...
                'OutputFormat','Samples');
            
            
%                 %Get randomized parameters
%                   Fc = randOverInterval(rangeFc);
%                 N = rangeNChip(randi(length(rangeNChip),1));
%                 Ncc = rangeNcc(randi(length(rangeNcc),1));

%                 %Get randomized parameters
                Fc = 1.75e+05;
                N = 4;
                Ncc = 1;
                
                % Create signal and update SNR
                chipWidth = Ncc/Fc;
                chipWidthSamples = round(chipWidth*fs)-1; % This must be an integer!
                chipWidth = chipWidthSamples*Ts;
                hPhaseZadoffChu.ChipWidth = chipWidth;
                hPhaseZadoffChu.NumChips = N;
                hPhaseZadoffChu.PRF = 1/((chipWidthSamples*N+1)*Ts);
                hPhaseZadoffChu.NumSamples = 1024;
                src = hPhaseZadoffChu();
                
% Pulse shape
% filter = phased.MatchedFilter( ...
%     'Coefficients',getMatchedFilter(hPhaseZadoffChu),...
%     'SampleRate', fs,...
%     'SpectrumWindow','None');
% tx = filter(src);  
filterCoeffs = rcosdesign(0.35, 4, sps);
tx = filter(filterCoeffs, 1, upsample(src,sps));

% Channel
channel = helperModClassTestChannel(...
  'SampleRate',fs, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / fs, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(tx);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(tx,rx)

% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction19,score19] = classify(trainedNet20SNR_1a,unknownFrames);

prediction19

plotScores(score19,modulationTypes)

%% BRPSKmodFrFtSignal
sampleFreq = 100e6;
sampleTime = 1/sampleFreq;
% sampleTime = 20.5e-9;
desSignalTime = 101.1e-6;
desSymFreq = 1e6;

rollOff = 0.5;
span = 2; % This must be an even number

redPhase = pi/12; % 20 degrees

f0 = 4e6;
f1 = 12e6;
chirpPhase = 0;
Ac = 1;

phi = -pi/4;
a = 2*phi/pi;

[sentSignal, timeVec, RRCtransSignal, freqSymbol] = BRPSKmodFrFTsignalGeneration(sampleTime, desSignalTime , desSymFreq, rollOff, span, redPhase, a);
% figure;
% plot(timeVec, RRCtransSignal);
% grid on;
% Channel
channel = helperModClassTestChannel(...
  'SampleRate',sampleFreq, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / sampleFreq, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(sentSignal);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(sentSignal,rx)

% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction20,score20] = classify(trainedNet20SNR_1a,unknownFrames);

prediction20

plotScores(score20,modulationTypes)

%% BRPSKmodChirpSignal

sampleFreq = 100e6;
sampleTime = 1/sampleFreq;
desChirpTime = 101.1e-6;
desSymFreq = 1e6;

rollOff = 0.5;
span = 2; % This must be an even number

redPhase = pi/12; % 20 degrees

f0 = 4e6;
f1 = 12e6;
chirpPhase = 0;
Ac = 1;

[sentSignal, timeVec, RRCtransSignal, complexChirp, freqSymbol] = BRPSKmodChirpSignalGeneration(sampleTime, desChirpTime, desSymFreq, rollOff, span, redPhase, f0, f1, chirpPhase, Ac);

% Square-root raised cosine filter
% filterCoeffs = rcosdesign(0.35, 4, 8);
% tx = filter(filterCoeffs, 1, upsample(sentSignal,8));

channel = helperModClassTestChannel(...
  'SampleRate',sampleFreq, ...
  'SNR',SNR, ...
  'PathDelays',[0 1.8 3.4] / sampleFreq, ...
  'AveragePathGains',[0 -2 -10], ...
  'KFactor',4, ...
  'MaximumDopplerShift',4, ...
  'MaximumClockOffset',5, ...
  'CenterFrequency',900e6);
rx = channel(RRCtransSignal);
% Plot transmitted and received signals
scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
  'LayoutDimensions',[2 1],'TimeSpan',45e-3);
scope(RRCtransSignal,rx)

% Frame generation for classification
unknownFrames = getNNFrames(rx,'Unknown');
% Classification
[prediction21,score21] = classify(trainedNet20SNR_1a,unknownFrames);

prediction21

plotScores(score21,modulationTypes)

%% Subroutines
function val = randOverInterval(interval)
% Expect interval to be <1x2> with format [minVal maxVal]
val = (interval(2) - interval(1)).*rand + interval(1);
end