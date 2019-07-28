function helperGenerateModulationfiles(wav,truth,Fs)
destinationFolder = 'C:\Users\racha\Documents\MATLAB\ModulationsDatabase';
 if ~exist(destinationFolder, 'dir')
  mkdir(destinationFolder);
end
%[~,~,~] = mkdir('ModulationsDatabase');
modTypes = unique(truth);

for idxM = 1:length(modTypes)
    modType = modTypes(idxM);
    [~,~,~] = mkdir(fullfile('ModulationsDatabase',char(modType)));
end

for idxW = 1:length(truth)
    sig = wav{idxW};
    if(isfinite(sig)==1)
      MOD = wvd(sig,Fs,'smoothedPseudo',kaiser(101,20),kaiser(101,20),'NumFrequencyPoints',500,'NumTimePoints',500);
      MOD = imresize(MOD,[227 227]);
      MOD = rescale(MOD);
      modType = truth(idxW);
      imwrite(MOD,fullfile('ModulationsDatabase',char(modType),sprintf('%s_%d.png',char(modType),idxW)))   
    end
    
    
end
end

