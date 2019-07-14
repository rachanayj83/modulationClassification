function plotSpectrogram(rxTest,rxTestLabel,modulationTypes,fs,sps)
%plotSpectrogram Spectrogram of frames

if size(rxTest,1) == 2
  IQAsRows = true;
else
  IQAsRows = false;
end
numRows = ceil(length(modulationTypes) / 4);
for modType=1:length(modulationTypes)
  subplot(numRows, 4, modType);
  idxOut = find(rxTestLabel == modulationTypes(modType), 1);
  if IQAsRows
    rxI = rxTest(1,:,1,idxOut);
    rxQ = rxTest(2,:,1,idxOut);
  else
    rxI = rxTest(1,:,1,idxOut);
    rxQ = rxTest(1,:,2,idxOut);
  end
  rx = squeeze(rxI) + 1i*squeeze(rxQ);
  spectrogram(rx,kaiser(sps),0,1024,fs,'centered');
  title(string(modulationTypes(modType)));
end
h = gcf; delete(findall(h.Children, 'Type', 'ColorBar'))
end
