function testAccuracy = sdrTest()
%sdrTest Test CNN performance with over-the-air signals
%   A = sdrTest sends test frames from one to another ADALM-PLUTO radio,
%   performs classification with the trained network, and returns the overall 
%   classification accuracy. Transmitting radio uses transmit-repeat capability
%   to send the same waveform repeatedly without loading the main loop.

modulationTypes = categorical(["BPSK", "QPSK", "8PSK", ...
  "2FSK","8FSK","GFSK", "CPFSK", ...
  "16QAM", "64QAM", "PAM4","LFM",...
  "Rect","Barker","Frank","P1","P2",...
  "P3","P4","Zadoff-Chu"]); 

load trainedNet20SNR_1a trainedNet
numFramesPerModType = 100;

sps = 8;                % Samples per symbol
spf = 1024;             % Samples per frame
fs = 1e6;             % Sample rate

txRadio = sdrtx('Pluto');
txRadio.RadioID = 'usb:0';
txRadio.CenterFrequency = 900e6;
txRadio.BasebandSampleRate = fs;

rxRadio = sdrrx('Pluto');
rxRadio.RadioID = 'usb:1';
rxRadio.CenterFrequency = 900e6;
rxRadio.BasebandSampleRate = fs;
rxRadio.SamplesPerFrame = spf;
rxRadio.ShowAdvancedProperties = true;
rxRadio.EnableQuadratureCorrection = false;
rxRadio.OutputDataType = 'single';

% Set random number generator to a known state to be able to regenerate
% the same frames every time the simulation is run
rng(1235)
tic

numModulationTypes = length(modulationTypes);
txModType = repmat(modulationTypes(1),numModulationTypes*numFramesPerModType,1);
estimatedModType = repmat(modulationTypes(1),numModulationTypes*numFramesPerModType,1);
frameCnt = 1;
for modType = 1:numModulationTypes
  fprintf('%s - Testing %s frames\n', ...
    datestr(toc/86400,'HH:MM:SS'), modulationTypes(modType))
  dataSrc = getSource(modulationTypes(modType), sps, 2*spf, fs);
  modulator = getModulator(modulationTypes(modType), sps, fs);
  if contains(char(modulationTypes(modType)), {'B-FM'})
    % Analog modulation types use a center frequency of 100 MHz
    txRadio.CenterFrequency = 100e6;
    rxRadio.CenterFrequency = 100e6;
    rxRadio.GainSource = 'Manual';
    rxRadio.Gain = 60;
  else
    % Digital modulation types use a center frequency of 900 MHz
    txRadio.CenterFrequency = 900e6;
    rxRadio.CenterFrequency = 900e6;
    rxRadio.GainSource = 'AGC Fast Attack';
  end
  
  % Start transmitter
  disp('Starting transmitter')
  x = dataSrc();
  y = modulator(x);
  y = y(4*sps+1:end,1);
  maxVal = max(max(abs(real(y))), max(abs(imag(y))));
  y = y *0.8/maxVal;
  % Download waveform signal to radio and repeatedly transmit it over the air
  transmitRepeat(txRadio, complex(y));
  
  disp('Starting receiver and test')
  for p=1:numFramesPerModType
    for frame=1:16
      rx = rxRadio();
    end
    
    frameEnergy = sum(abs(rx).^2);
    rx = rx / sqrt(frameEnergy);
    reshapedRx(1,:,1,1) = real(rx);
    reshapedRx(1,:,2,1) = imag(rx);
    
    % Classify
    txModType(frameCnt) = modulationTypes(modType);
    estimatedModType(frameCnt) = classify(trainedNet, reshapedRx);
    
    frameCnt = frameCnt + 1;
    
    % Pause for 0.1 seconds to get an independent channel (assuming a
    % channel coherence time of less than 0.1 seconds)
    pause(0.1)
  end
  disp('Releasing radios')
  release(txRadio);
  release(rxRadio);
end
testAccuracy = mean(txModType == estimatedModType);
disp("Test accuracy: " + testAccuracy*100 + "%")

figure
cm = confusionchart(txModType, estimatedModType);
cm.Title = 'Confusion Matrix for Test Data';
cm.RowSummary = 'row-normalized';
cm.Parent.Position = [cm.Parent.Position(1:2) 740 424];
end

