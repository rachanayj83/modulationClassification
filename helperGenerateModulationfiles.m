function helperGenerateModulationfiles(wav,truth,Fs)
[~,~,~] = mkdir('ModulationsDatabase');
modTypes = unique(truth);

for idxM = 1:length(modTypes)
    modType = modTypes(idxM);
    [~,~,~] = mkdir(fullfile('ModulationsDatabase',char(modType)));
end

for idxW = 1:length(truth)
    sig = wav{idxW};
    MOD = wvd(sig,Fs,'smoothedPseudo',kaiser(101,20),kaiser(101,20),'NumFrequencyPoints',500,'NumTimePoints',500);
    MOD = imresize(MOD,[227 227]);
    MOD = rescale(MOD);
    modType = truth(idxW);
    
    imwrite(MOD,fullfile('ModulationsDatabase',char(modType),sprintf('%d.png',idxW)))
end
end
