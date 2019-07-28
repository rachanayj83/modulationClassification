%clc
%close all
%restoredefaultpath %restores the default path of the functions
clear all;

%% 
%addpath('fracft/'); %The path were the functions are
%addpath('Functions/');

% The sample freq must be bigger than twice the biggest frequency 
% of the chirp signal, because of the Nyquist sampling theorem, aliasing
% and to have a good estimation of the SCF. See what happens when you
% change from 30MHz to 25MHz. In principle, 25 should be enough to have a
% good estimation of the SCF, as you can see up to 12.5MHz, but it gives a
% poor SCF. Whereas, 30MHz, with which you can see up to 15MHz, gives a
% good SCF. Also to see well the cyclic freq your sampling freq should be a
% multiple of the symbol freq.

sampleFreq = 100e6;
sampleTime = 1/sampleFreq;
% sampleTime = 20.5e-9;
desSignalTime = 101.1e-6;
desSymFreq = 1e6;

rollOff = 0.5;
span = 2; % This must be an even number

redPhase = pi/12; % 20 degrees

f0 = 4e6;
f1 = 12e6;
chirpPhase = 0;
Ac = 1;

phi = -pi/4;
a = 2*phi/pi;

[sentSignal, timeVec, RRCtransSignal, freqSymbol] = BRPSKmodFrFTsignalGeneration(sampleTime, desSignalTime , desSymFreq, rollOff, span, redPhase, a);

figure;
plot(timeVec, sentSignal);
grid on;

%% SCF
alphaLimit = 10e6;
alphaResol = 50e3;
freqSCFresol = 250e3;
restrictFactor = 39/40;

[SCF, ~, ~, freq, alphaFreq, ~] = RestrictedSCFcomplex(sentSignal, alphaResol, alphaLimit, freqSCFresol, sampleTime, restrictFactor);

% figure()
% s = surf(freq/1e6, alphaFreq/1e6, (abs(SCF)));
% s.EdgeColor = 'none';
% xlabel('f (MHz)')
% ylabel('\alpha (MHz)')

%% Projection of the alpha frequencies of the SCF
% figure;
% plot(alphaFreq, db(mean(abs((SCF)), 2)));
% grid on;