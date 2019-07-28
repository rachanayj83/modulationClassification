[wav, modType] = helperRadarWaveformGeneration();

idLFM = find(modType == "LFM",3);
nfft = 2^nextpow2(length(wav{1}));
f = (0:(nfft/2-1))/nfft*100e6;

figure
subplot(1,3,1)
Z = fft(wav{idLFM(1)},nfft);
plot(f/1e6,abs(Z(1:nfft/2)))
xlabel('Frequency (MHz)');ylabel('Amplitude');axis square
subplot(1,3,2)
Z = fft(wav{idLFM(2)},nfft);
plot(f/1e6,abs(Z(1:nfft/2)))
xlabel('Frequency (MHz)');ylabel('Amplitude');axis square
subplot(1,3,3)
Z = fft(wav{idLFM(3)},nfft);
plot(f/1e6,abs(Z(1:nfft/2)))
xlabel('Frequency (MHz)');ylabel('Amplitude');axis square
%% 
figure
subplot(3,3,1)
wvd(wav{find(modType == "Rect",1)},100e6,'smoothedPseudo')
axis square; colorbar off; title('Rect')
subplot(3,3,2)
wvd(wav{find(modType == "LFM",1)},100e6,'smoothedPseudo')
axis square; colorbar off; title('LFM')
subplot(3,3,3)
wvd(wav{find(modType == "Barker",1)},100e6,'smoothedPseudo')
axis square; colorbar off; title('Barker')

subplot(3,3,4)
wvd(wav{find(modType == "Frank",1)},100e6,'smoothedPseudo')
axis square; colorbar off; title('Frank')
subplot(3,3,5)
wvd(wav{find(modType == "P1",1)},100e6,'smoothedPseudo')
axis square; colorbar off; title('P1')
subplot(3,3,6)
wvd(wav{find(modType == "P2",1)},100e6,'smoothedPseudo')
axis square; colorbar off; title('P2')

subplot(3,3,7)
wvd(wav{find(modType == "P3",1)},100e6,'smoothedPseudo')
axis square; colorbar off; title('P3')
subplot(3,3,8)
wvd(wav{find(modType == "P4",1)},100e6,'smoothedPseudo')
axis square; colorbar off; title('P4')
% subplot(3,3,9)
% wvd(wav{find(modType == "ZadoffChu",1)},100e6,'smoothedPseudo')
% axis square; colorbar off; title('ZadoffChu')


%% 
 helperGenerateModulationfiles(wav,modType,100e6)
 
 %%
folders = fullfile('ModulationsDatabase',{'Rect','LFM','Barker','Frank','P1','P2','P3','P4','ZadoffChu'});