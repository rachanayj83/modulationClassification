function modulator = getModulator(modType, sps, fs)
%getModulator Modulation function selector
%   MOD = getModulator(TYPE,SPS,FS) returns the modulator function handle
%   MOD based on TYPE. SPS is the number of samples per symbol and FS is 
%   the sample rate.

switch modType
  case "BPSK"
    modulator = @(x)bpskModulator(x,sps);
  case "QPSK"
    modulator = @(x)qpskModulator(x,sps);
  case "8PSK"
    modulator = @(x)psk8Modulator(x,sps);
  case "2FSK"
    modulator = @(x)fsk2Modulator(x,sps);
  case "8FSK"
    modulator = @(x)fsk8Modulator(x,sps);
  case "GFSK"
    modulator = @(x)gfskModulator(x,sps);
  case "CPFSK"
    modulator = @(x)cpfskModulator(x,sps);
  case "16QAM"
    modulator = @(x)qam16Modulator(x,sps);
  case "64QAM"
    modulator = @(x)qam64Modulator(x,sps);
  case "PAM4"
    modulator = @(x)pam4Modulator(x,sps);
  case "B-FM"
    modulator = @(x)bfmModulator(x, fs);
%   case "DSB-AM"
%     modulator = @(x)dsbamModulator(x, fs);
%   case "SSB-AM"
%     modulator = @(x)ssbamModulator(x, fs);
   case "LFM"
     modulator = @(x)lfmModulator(x);
   case "Rect"
     modulator = @(x)rectModulator(x);
   case "Barker"
     modulator = @(x)barkerModulator(x);
   case "Frank"
    modulator = @(x)frankModulator(x);
   case "P1"
    modulator = @(x)p1Modulator(x);
   case "P2"
    modulator = @(x)p2Modulator(x);
   case "P3"
    modulator = @(x)p3Modulator(x);
   case "P4"
    modulator = @(x)p4Modulator(x);
   case "Zadoff-Chu"
    modulator = @(x)zadoffchuModulator(x);
end

%% Modulators

function y = bpskModulator(x,sps)
%bpskModulator BPSK modulator with pulse shaping
%   Y=bpskModulator(X,SPS) BPSK modulates input, X, and returns a root-raised
%   cosine pulse shaped signal, Y. X must be a column vector of values in
%   the set of [0 1]. Root-raised cosine filter has a roll-off factor of
%   0.35 and spans four symbols. The output signal, Y, has unit power.

persistent filterCoeffs
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
end
% Modulate
syms = pskmod(x,2);
% Pulse shape
y = filter(filterCoeffs, 1, upsample(syms,sps));
end

function y = qpskModulator(x,sps)
%qpskModulator QPSK modulator with pulse shaping
%   Y=qpskModulator(X,SPS) QPSK modulates input, X, and returns a root-raised
%   cosine pulse shaped signal, Y. X must be a column vector of values in
%   the set of [0 3]. Root-raised cosine filter has a roll-off factor of
%   0.35 and spans four symbols. The output signal, Y, has unit power.

persistent filterCoeffs
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
end
% Modulate
syms = pskmod(x,4,pi/4);
% Pulse shape
y = filter(filterCoeffs, 1, upsample(syms,sps));
end

function y = psk8Modulator(x,sps)
%psk8Modulator 8-PSK modulator with pulse shaping
%   Y = psk8Modulator(X,SPS) 8-PSK modulates the input X, and returns the 
%   root-raised cosine pulse shaped signal Y. X must be a column vector 
%   of values in the set [0 7]. The root-raised cosine filter has a 
%   roll-off factor of 0.35 and spans four symbols. The output signal 
%   Y has unit power.

persistent filterCoeffs
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
end
% Modulate
syms = pskmod(x,8);
% Pulse shape
y = filter(filterCoeffs, 1, upsample(syms,sps));
end

function y = fsk2Modulator(x,sps)
%fsk2Modulator FSK modulator with pulse shaping
%   Y = fsk2Modulator(X,SPS) FSK modulates the input X, and returns the 
%   root-raised cosine pulse shaped signal Y. X must be a column vector 
%   of values in the set [0 1]. The root-raised cosine filter has a 
%   roll-off factor of 0.35 and spans four symbols. The output signal 
%   Y has unit power.

persistent mod 
if isempty(mod)
  M = 2;        % Modulation order
  freqsep = 4;  % Frequency separation (Hz)
  symbolRate = 45;
  %fsamp = sps*symbolRate;
  
 mod = comm.FSKModulator(...
    'ModulationOrder', M, ...
    'SymbolMapping', 'Binary', ...
    'FrequencySeparation', freqsep, ...
    'SymbolRate', symbolRate, ...
    'SamplesPerSymbol', sps);
end
% Modulate
y = mod(x);
end


function y = fsk8Modulator(x,sps)
%fsk8Modulator FSK modulator with pulse shaping
%   Y = fsk8Modulator(X,SPS) FSK modulates the input X, and returns the 
%   root-raised cosine pulse shaped signal Y. X must be a column vector 
%   of values in the set [0 1]. The root-raised cosine filter has a 
%   roll-off factor of 0.35 and spans four symbols. The output signal 
%   Y has unit power.

persistent mod 
if isempty(mod)
  M = 8;        % Modulation order
  freqsep = 4;  % Frequency separation (Hz)
  symbolRate = 45;
  
  mod = comm.FSKModulator(...
    'ModulationOrder', M, ...
    'FrequencySeparation', freqsep, ...
    'SymbolRate', symbolRate, ...
    'SamplesPerSymbol', sps);
end 
% Modulate
y = mod(x);
end


function y = gfskModulator(x,sps)
%gfskModulator GFSK modulator
%   Y=gfskModulator(X,SPS) GFSK modulates input, X, and returns signal, Y. 
%   X must be a column vector of values in the set of [0 1]. BT product is
%   0.35 and modulation index is 1. The output signal, Y, has unit power.

persistent mod meanM
if isempty(mod)
  M = 2;
  mod = comm.CPMModulator(...
    'ModulationOrder', M, ...
    'FrequencyPulse', 'Gaussian', ...
    'BandwidthTimeProduct', 0.35, ...
    'ModulationIndex', 1, ...
    'SamplesPerSymbol', sps);
  meanM = mean(0:M-1);
end
% Modulate
y = mod(2*(x-meanM));
end

function y = cpfskModulator(x,sps)
%cpfskModulator CPFSK modulator
%   Y=cpfskModulator(X,SPS) CPFSK modulates input, X, and returns signal, Y. 
%   X must be a column vector of values in the set of [0 1]. Modulation 
%   index is 0.5. The output signal, Y, has unit power.

persistent mod meanM
if isempty(mod)
  M = 2;
  mod = comm.CPFSKModulator(...
    'ModulationOrder', M, ...
    'ModulationIndex', 0.5, ...
    'SamplesPerSymbol', sps);
  meanM = mean(0:M-1);
end
% Modulate
y = mod(2*(x-meanM));
end

function y = qam16Modulator(x,sps)
%qam16Modulator 16-QAM modulator with pulse shaping
%   Y = qam16Modulator(X,SPS) 16-QAM modulates the input X, and returns the 
%   root-raised cosine pulse shaped signal Y. X must be a column vector 
%   of values in the set [0 15]. The root-raised cosine filter has a 
%   roll-off factor of 0.35 and spans four symbols. The output signal 
%   Y has unit power.

persistent filterCoeffs
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
end
% Modulate and pulse shape
syms = qammod(x,16,'UnitAveragePower',true);
% Pulse shape
y = filter(filterCoeffs, 1, upsample(syms,sps));
end

function y = qam64Modulator(x,sps)
%qam64Modulator 64-QAM modulator with pulse shaping
%   Y = qam64Modulator(X,SPS) 64-QAM modulates the input X, and returns the 
%   root-raised cosine pulse shaped signal Y. X must be a column vector 
%   of values in the set [0 63]. The root-raised cosine filter has a 
%   roll-off factor of 0.35 and spans four symbols. The output signal 
%   Y has unit power.

persistent filterCoeffs
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
end
% Modulate
syms = qammod(x,64,'UnitAveragePower',true);
% Pulse shape
y = filter(filterCoeffs, 1, upsample(syms,sps));
end


function y = pam4Modulator(x,sps)
%pam4Modulator PAM4 modulator with pulse shaping
%   Y=qam16Modulator(X,SPS) PAM4 modulates input, X, and returns a root-raised
%   cosine pulse shaped signal, Y. X must be a column vector of values in
%   the set of [0 3]. Root-raised cosine filter has a roll-off factor of
%   0.35 and spans four symbols. The output signal, Y, has unit power.

persistent filterCoeffs amp
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
  amp = 1 / sqrt(mean(abs(pammod(0:3, 4)).^2));
end
% Modulate
syms = amp * pammod(x,4);
% Pulse shape
y = filter(filterCoeffs, 1, upsample(syms,sps));
end

function y = lfmModulator(x)
persistent filterCoeffs
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
end

y = filter(filterCoeffs, 1,upsample(x,sps));
end

function y = rectModulator(x)
 persistent filterCoeffs
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
end
y = filter(filterCoeffs, 1,upsample(x,sps));
end

function y = barkerModulator(x)
  persistent filterCoeffs
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
end
y = filter(filterCoeffs, 1,upsample(x,sps));
end

function y = frankModulator(x)
persistent filterCoeffs
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
end
y = filter(filterCoeffs, 1,upsample(x,sps));
end

function y = p1Modulator(x)
   persistent filterCoeffs
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
end
y = filter(filterCoeffs, 1,upsample(x,sps));
end

function y = p2Modulator(x)
persistent filterCoeffs
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
end
y = filter(filterCoeffs, 1, upsample(x,sps));
end

function y = p3Modulator(x)
persistent filterCoeffs
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
end
y = filter(filterCoeffs, 1, upsample(x,sps));
end

function y = p4Modulator(x)
persistent filterCoeffs
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
end
y = filter(filterCoeffs, 1, upsample(x,sps));
end

function y = zadoffchuModulator(x)
persistent filterCoeffs
if isempty(filterCoeffs)
  filterCoeffs = rcosdesign(0.35, 4, sps);
end
y = filter(filterCoeffs, 1, upsample(x,sps));
end


%% Subroutines
function val = randOverInterval(interval)
% Expect interval to be <1x2> with format [minVal maxVal]
val = (interval(2) - interval(1)).*rand + interval(1);
end
end




