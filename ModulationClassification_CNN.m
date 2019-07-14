% clear all;
% clc;
%%
modulationTypes = categorical(["BPSK", "QPSK", "8PSK", ...
  "2FSK","8FSK","GFSK", "CPFSK", ...
  "16QAM", "64QAM", "PAM4","LFM",...
  "Rect","Barker","Frank","P1","P2",...
  "P3","P4","Zadoff-Chu"]); 

%   "B-FM", "DSB-AM", "SSB-AM"]);

%% Generating Waveforms for training

numFramesPerModType = 15000;
percentTrainingSamples = 80;
percentValidationSamples = 10;
percentTestSamples = 10;

sps = 8;                % Samples per symbol
spf = 1024;             % Samples per frame
symbolsPerFrame = spf / sps;
fs = 1e6;             % Sample rate
% fs = 200e3;             % Sample rate
fc = [900e6 100e6];     % Center frequencies

% Initilisations for RADAR waveforms
rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
% rangeN = [512, 1920]; % Number of collected signal samples range
% rangeB = [fs/20, fs/16]; % Bandwidth (Hz) range
% sweepDirections = {'Up','Down'}






%% Create Channel Impairments
% Pass each frame through a channel with AWGN,Rician multipath fading,
% Clock offset, resulting in center frequency offset and sampling time drift

SNR = -5;
std = sqrt(10.^(-SNR/-5))

awgnChannel = comm.AWGNChannel(...
  'NoiseMethod', 'Signal to noise ratio (SNR)', ...
  'SignalPower', 1, ...
  'SNR', SNR)

%% Rician Multipath
% The channel passes the signals through a Rician multipath fading channel using the 
% comm.RicianChannel System object. Assume a delay profile of [0 1.8 3.4] samples with 
% corresponding average path gains of [0 -2 -10] dB. The K-factor is 4 and the maximum
% Doppler shift is 4 Hz

multipathChannel = comm.RicianChannel(...
  'SampleRate', fs, ...
  'PathDelays', [0 1.8 3.4]/fs, ...
  'AveragePathGains', [0 -2 -10], ...
  'KFactor', 4, ...
  'MaximumDopplerShift', 4)

%% Clock Offset

maxDeltaOff = 5;
deltaOff = (rand()*2*maxDeltaOff) - maxDeltaOff;
C = 1 + (deltaOff/1e6);

%% Frequency Offset

% Subject each frame to a frequency offset based on clock offset factor C and the center 
% frequency

offset = -(C-1)*fc(1);
frequencyShifter = comm.PhaseFrequencyOffset(...
  'SampleRate', fs, ...
  'FrequencyOffset', offset)

%% Sampling Rate Offset

%Subject each frame to a sampling rate offset based on clock offset factor C

channel = helperModClassTestChannel(...
  'SampleRate', fs, ...
  'SNR', SNR, ...
  'PathDelays', [0 1.8 3.4] / fs, ...
  'AveragePathGains', [0 -2 -10], ...
  'KFactor', 4, ...
  'MaximumDopplerShift', 4, ...
  'MaximumClockOffset', 5, ...
  'CenterFrequency', 900e6)

chInfo = info(channel)


%% Waveform Generation

% Set the random number generator to a known state to be able to regenerate
% the same frames every time the simulation is run
rng(1235)
tic

numModulationTypes = length(modulationTypes);

channelInfo = info(channel);
frameStore = helperModClassFrameStore(...
  numFramesPerModType*numModulationTypes,spf,modulationTypes);
transDelay = 50;
for modType = 1:numModulationTypes
  fprintf('%s - Generating %s frames\n', ...
    datestr(toc/86400,'HH:MM:SS'), modulationTypes(modType))
  numSymbols = (numFramesPerModType / sps);
  dataSrc = getSource(modulationTypes(modType), sps, 2*spf, fs);
  modulator = getModulator(modulationTypes(modType), sps, fs, SNR);
  if contains(char(modulationTypes(modType)), {'B-FM','DSB-AM','SSB-AM'})
    % Analog modulation types use a center frequency of 100 MHz
    channel.CenterFrequency = 100e6;
  else
    % Digital modulation types use a center frequency of 900 MHz
    channel.CenterFrequency = 900e6;
  end
  
  for p=1:numFramesPerModType
    % Generate random data
    x = dataSrc();
    
    % Modulate
    y = modulator(x);
    
    if contains(char(modulationTypes(modType)), {'LFM'})
      release(channel);
    end
    % Pass through independent channels
    rxSamples = channel(y);
    
    % Remove transients from the beginning, trim to size, and normalize
    frame = helperModClassFrameGenerator(rxSamples, spf, spf, transDelay, sps);
    
    % Add to frame store
    add(frameStore, frame, modulationTypes(modType));
  end
end

%% 
% Next divide the frames into training, validation, and test data. 
% By default, frameStore places I/Q baseband samples in rows in the output frames. 
% The output frames have the size [2xspfx1xN], where the first row is in-phase samples 
% and the second row is quadrature samples.

[mcfsTraining,mcfsValidation,mcfsTest] = splitData(frameStore,...
  [percentTrainingSamples,percentValidationSamples,percentTestSamples]);
% [rxTraining,rxTrainingLabel] = get(mcfsTraining);
% [rxValidation,rxValidationLabel] = get(mcfsValidation);
% [rxTest,rxTestLabel] = get(mcfsTest);

% Put the data in [1xspfx2] format
mcfsTraining.OutputFormat = "IQAsPages";
[rxTraining,rxTrainingLabel] = get(mcfsTraining);
mcfsValidation.OutputFormat = "IQAsPages";
[rxValidation,rxValidationLabel] = get(mcfsValidation);
mcfsTest.OutputFormat = "IQAsPages";
[rxTest,rxTestLabel] = get(mcfsTest);

%% 

% Plot the amplitude of the real and imaginary parts of the example frames
% against the sample number
figure
plotTimeDomain(rxTest,rxTestLabel,modulationTypes,fs)

%%
% Plot a spectrogram of the example frames
% figure
% plotSpectrogram(rxTest,rxTestLabel,modulationTypes,fs,sps)

%%
% Plot the label distributions
% figure
% subplot(3,1,1)
% histogram(rxTrainingLabel)
% title("Training Label Distribution")
% subplot(3,1,2)
% histogram(rxValidationLabel)
% title("Validation Label Distribution")
% subplot(3,1,3)
% histogram(rxTestLabel)
% title("Test Label Distribution")

%% Train the CNN

dropoutRate = 0.6;
numModTypes = numel(modulationTypes);
netWidth = 1;
filterSize = [1 sps];
poolSize = [1 2];
 % Set the input layer input size to [1xspfx2]
modClassNet = [
  imageInputLayer([1 spf 2], 'Normalization', 'none', 'Name', 'Input Layer')
  
  convolution2dLayer(filterSize, 16*netWidth, 'Padding', 'same', 'Name', 'CNN1')
  batchNormalizationLayer('Name', 'BN1')
  reluLayer('Name', 'ReLU1')
  maxPooling2dLayer(poolSize, 'Stride', [1 2], 'Name', 'MaxPool1')
  
  convolution2dLayer(filterSize, 24*netWidth, 'Padding', 'same', 'Name', 'CNN2')
  batchNormalizationLayer('Name', 'BN2')
  reluLayer('Name', 'ReLU2')
  maxPooling2dLayer(poolSize, 'Stride', [1 2], 'Name', 'MaxPool2')
  
  convolution2dLayer(filterSize, 32*netWidth, 'Padding', 'same', 'Name', 'CNN3')
  batchNormalizationLayer('Name', 'BN3')
  reluLayer('Name', 'ReLU3')
  maxPooling2dLayer(poolSize, 'Stride', [1 2], 'Name', 'MaxPool3')
  
  convolution2dLayer(filterSize, 48*netWidth, 'Padding', 'same', 'Name', 'CNN4')
  batchNormalizationLayer('Name', 'BN4')
  reluLayer('Name', 'ReLU4')
  maxPooling2dLayer(poolSize, 'Stride', [1 2], 'Name', 'MaxPool4')
  
  convolution2dLayer(filterSize, 64*netWidth, 'Padding', 'same', 'Name', 'CNN5')
  batchNormalizationLayer('Name', 'BN5')
  reluLayer('Name', 'ReLU5')
  maxPooling2dLayer(poolSize, 'Stride', [1 2], 'Name', 'MaxPool5')
  
  convolution2dLayer(filterSize, 96*netWidth, 'Padding', 'same', 'Name', 'CNN6')
  batchNormalizationLayer('Name', 'BN6')
  reluLayer('Name', 'ReLU6')
  
  convolution2dLayer(filterSize, 128*netWidth, 'Padding', 'same', 'Name', 'CNN7')
  batchNormalizationLayer('Name', 'BN7')
  reluLayer('Name', 'ReLU7')
  
  averagePooling2dLayer([1 ceil(spf/32)], 'Name', 'AP1')
  
  fullyConnectedLayer(numModTypes, 'Name', 'FC1')
  softmaxLayer('Name', 'SoftMax')
  
  classificationLayer('Name', 'Output') ]

analyzeNetwork(modClassNet)

%% Setting training options 

maxEpochs = 15;
miniBatchSize = 256;
validationFrequency = floor(numel(rxTrainingLabel)/miniBatchSize);
options = trainingOptions('adam', ...
  'InitialLearnRate',1e-2, ...
  'MaxEpochs',maxEpochs, ...
  'MiniBatchSize',miniBatchSize, ...
  'Shuffle','every-epoch', ...
  'Plots','training-progress', ...
  'Verbose',false, ...
  'ValidationData',{rxValidation,rxValidationLabel}, ...
  'ValidationFrequency',validationFrequency, ...
  'LearnRateSchedule', 'piecewise', ...
  'LearnRateDropPeriod', 9, ...
  'LearnRateDropFactor', 0.1, ...
  'ExecutionEnvironment', 'multi-gpu');
%%

trainNow = true;
if trainNow == true
  fprintf('%s - Training the network\n', datestr(toc/86400,'HH:MM:SS'))
  trainedNet0SNR_v7 = trainNetwork(rxTraining,rxTrainingLabel,modClassNet,options);
else
  load trainedNet20SNR_v1
end

%% Evaluate Trained network

fprintf('%s - Classifying test frames\n', datestr(toc/86400,'HH:MM:SS'))
rxTestPred = classify(trainedNet0SNR_v7,rxTest);
testAccuracy = mean(rxTestPred == rxTestLabel);
disp("Test accuracy: " + testAccuracy*100 + "%")

%% Plot confusion matrix

figure
cm = confusionchart(rxTestLabel, rxTestPred);
cm.Title = 'Confusion Matrix for Test Data';
cm.RowSummary = 'row-normalized';
% cm.Normalization = 'total-normalized';
sortClasses(cm,'descending-diagonal')
cm.Parent.Position = [cm.Parent.Position(1:2) 740 424];


