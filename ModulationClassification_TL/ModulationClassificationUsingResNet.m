 [wav, modType] = helperRadarWaveformGeneration();
 
 %% 
 helperGenerateModulationfiles(wav,modType,100e6)
 
 %%
 folders = fullfile('ModulationsDatabase',{'Rect','LFM','Barker','Frank','P1','P2','P3','P4','Zadoff-Chu'});

%%
[wav, modType] = helperCommsWaveformGeneration();

%% 
 helperGenerateModulationfiles(wav,modType,100e6)

%%
folders = fullfile('ModulationsDatabase',{'Rect','LFM','Barker','Frank','P1','P2','P3','P4','Zadoff-Chu','GFSK','CPFSK','BPSK','2FSK','8FSK','QPSK','8PSK','16QAM','64QAM','PAM4'});

imds = imageDatastore(folders,...
    'FileExtensions','.png','LabelSource','foldernames','ReadFcn',@readModulationsDatabaseForResNet);
%% Splitting the data for training, testing snd validation
[imdsTrain,imdsTest,imdsValidation] = splitEachLabel(imds,0.8,0.1);

%% Define network architecture

numClasses = 19;
net = resnet50;
layersTransfer = net.Layers(1:end-3);
layers = [
    layersTransfer
    fullyConnectedLayer(numClasses,'WeightLearnRateFactor',20,'BiasLearnRateFactor',20)
    softmaxLayer
    classificationLayer]

options = trainingOptions('sgdm', ...
    'MiniBatchSize',256, ...
    'MaxEpochs',5, ...
    'InitialLearnRate',1e-3, ...
    'Shuffle','every-epoch', ...
    'Verbose',false, ...
    'Plots','training-progress',...
    'ValidationData',imdsValidation);

%% Train the network

trainedNet = trainNetwork(shuffle(imdsTrain),layers,options);

%% Evaluate Performance on the signals

predicted = classify(trainedNet,imdsTest);
figure;
confusionchart(predicted,imdsTest.Labels,'Normalization','column-normalized');

% %%
% GFSK_GFSK = readimage(imdsTest,find((imdsTest.Labels == 'GFSK') & (predicted == 'GFSK'),1));
% GFSK_2FSK = readimage(imdsTest,find((imdsTest.Labels == 'GFSK') & (predicted == '2FSK'),1));
% PAM4_PAM4 = readimage(imdsTest,find((imdsTest.Labels == 'PAM4') & (predicted == 'PAM4'),1));
% PAM4_BPSK = readimage(imdsTest,find((imdsTest.Labels == 'PAM4') & (predicted == 'BPSK'),1));
% 
% figure
% subplot(2,2,1)
% imagesc(GFSK_GFSK(:,:,1))
% axis square; title({'Actual Class: GFSK','Predicted Class: GFSK'})
% subplot(2,2,2)
% imagesc(GFSK_2FSK(:,:,1))
% axis square; title({'Actual Class: GFSK','Predicted Class: 2FSK'})
% subplot(2,2,3)
% imagesc(PAM4_PAM4(:,:,1))
% axis square; title({'Actual Class: PAM4','Predicted Class: PAM4'})
% subplot(2,2,4)
% imagesc(PAM4_BPSK(:,:,1))
% axis square; title({'Actual Class: PAM4','Predicted Class: BPSK'})
% 
% %% Cross-validate trained network
% load trainednet
% 
% %%
% % Set the random number generator to a known state to be able to regenerate
% % the same frames every time the simulation is run
% rng(123456)
% fs = 1e8;
% rangeFc = [fs/6, fs/5]; % Center frequency (Hz) range
% % Random bits
% d = randi([0 1],64,1);
% % BPSK modulation
% syms = pskmod(d,2);
% % Square-root raised cosine filter
% filterCoeffs = rcosdesign(0.35,4,8);
% tx = filter(filterCoeffs,1,upsample(syms,8));
% % Channel
% multipathChannel = comm.RicianChannel(...
%     'SampleRate', fs, ...
%     'PathDelays', [0 1.8 3.4]/fs, ...
%     'AveragePathGains', [0 -2 -10], ...
%     'KFactor', 4, ...
%     'MaximumDopplerShift', 4);
% 
% hFreqOffset = comm.PhaseFrequencyOffset(...
%     'SampleRate',fs);
% 
% rx = multipathChannel(tx);
% 
% wav = awgn(rx,20);
%                 
%                  % Add frequency offset
%                  fc = randOverInterval(rangeFc);
%                  hFreqOffset.FrequencyOffset = 900e3;
%                  wav = hFreqOffset(wav); % Frequency shift
%                  
%                  rx = multipathChannel(wav);
%                  
% % Plot transmitted and received signals
% scope = dsp.TimeScope(2,200e3,'YLimits',[-1 1],'ShowGrid',true,...
%   'LayoutDimensions',[2 1],'TimeSpan',45e-3);
% scope(tx,rx)
% 
% %% Cross-validate trained network
% MOD_Predict = wvd(rx,fs,'smoothedPseudo',kaiser(101,20),kaiser(101,20),'NumFrequencyPoints',500,'NumTimePoints',500);
% 
% MOD_Predict = imresize(MOD_Predict,[227 227]);
% MOD_Predict = rescale(MOD_Predict);
% MOD_Predict = repmat(MOD_Predict, 1, 1, 3);  
% % imwrite(MOD_Predict,fullfile('ModulationsPredictDatabase','Unknown',sprintf('Unknown.png')));
% % Classification
% predicted_modulation = char(classify(trainedNet,MOD_Predict));
% 
% %%
% 
% predicted_modulation

%% Subroutines
function val = randOverInterval(interval)
% Expect interval to be <1x2> with format [minVal maxVal]
val = (interval(2) - interval(1)).*rand + interval(1);
end