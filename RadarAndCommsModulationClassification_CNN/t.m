function plotTimeDomain(rxTest,rxTestLabel,modulationTypes,fs)
%plotTimeDomain Time domain plots of frames

numRows = ceil(length(modulationTypes) / 4);
spf = size(rxTest,2);
t = 1000*(0:spf-1)/fs;
if size(rxTest,1) == 2
  IQAsRows = true;
else
  IQAsRows = false;
end
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
  plot(t,squeeze(rxI), '-'); grid on; axis equal; axis square
  hold on
  plot(t,squeeze(rxQ), '-'); grid on; axis equal; axis square
  hold off
  title(string(modulationTypes(modType)));
  xlabel('Time (ms)'); ylabel('Amplitude')
end
end