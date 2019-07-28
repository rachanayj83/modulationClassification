clear all;
clc;

modulationTypes = categorical(["LFM","Rect","Barker","Frank","P1","P2","P3","P4","Zadoff-Chu",...
    "BPSK","QPSK","8PSK","PAM4","16QAM","64QAM",...
    "GFSK","CPFSK","2FSK","8FSK"]);
    
numModulationTypes = length(modulationTypes);
numFramesPerModType = 10;

percentTrainingSamples = 80;
percentValidationSamples = 10;
percentTestSamples = 10;

sps = 8;                % Samples per symbol
spf = 1024;             % Samples per frame
symbolsPerFrame = spf / sps;
fs = 100e3;             % Sample rate
fc = [900e6 100e6];     % Center frequencies



frameStore = helperModClassFrameStore(...
  numFramesPerModType*numModulationTypes,spf,modulationTypes);

helperModClassification_CNN(frameStore);

%%
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
 figure
 plotSpectrogram(rxTest,rxTestLabel,modulationTypes,fs,sps)
%%
% Plot the label distributions
figure
subplot(3,1,1)
histogram(rxTrainingLabel)
title("Training Label Distribution")
subplot(3,1,2)
histogram(rxValidationLabel)
title("Validation Label Distribution")
subplot(3,1,3)
histogram(rxTestLabel)
title("Test Label Distribution")

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
  trainedNet20SNR_v1a = trainNetwork(rxTraining,rxTrainingLabel,modClassNet,options);
else
  load trainedNet20SNR_v1
end

%% Evaluate Trained network

fprintf('%s - Classifying test frames\n', datestr(toc/86400,'HH:MM:SS'))
rxTestPred = classify(trainedNet20SNR_v1a,rxTest);
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