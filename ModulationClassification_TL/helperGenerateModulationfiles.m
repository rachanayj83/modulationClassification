function helperGenerateModulationfiles(wav,truth,Fs)
if ~exist('ModulationsDatabase', 'dir')
 [~,~,~] = mkdir('ModulationsDatabase');
end

modTypes = unique(truth);

for idxM = 1:length(modTypes)
    modType = modTypes(idxM);
    [~,~,~] = mkdir(fullfile('ModulationsDatabase',char(modType)));
end

for idxW = 1:length(truth)
    sig = wav{idxW};
    MOD = wvd(sig,Fs,'smoothedPseudo',kaiser(101,20),kaiser(101,20),'NumFrequencyPoints',500,'NumTimePoints',500);
%     MOD = ind2rgb(im2uint8(rescale(MOD)),jet(128));
%For AlexNet image size is [227 227]
%     MOD = imresize(MOD,[227 227]);
    MOD = imresize(MOD,[224 224]);
    MOD = rescale(MOD);
    modType = truth(idxW);
    
    imwrite(MOD,fullfile('ModulationsDatabase',char(modType),sprintf('%s_%d.png',char(modType),idxW)))
end
end
