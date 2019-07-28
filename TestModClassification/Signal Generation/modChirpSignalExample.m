%clc
%close all
%restoredefaultpath %restores the default path of the functions
clear all;

%% 
addpath('Functions/'); %The path were the functions are

% The sample freq must be bigger than twice the biggest frequency 
% of the chirp signal, because of the Nyquist sampling theorem, aliasing
% and to have a good estimation of the SCF. See what happens when you
% change from 46MHz to 25MHz. In principle, 25 should be enough to have a
% good estimation of the SCF, as you can see up to 12.5MHz, but it gives a
% poor SCF. Whereas, 46MHz, with which you can see up to 23MHz, gives a
% good SCF. Also to be able to see well the cyclic freq your sampling freq 
% should be a multiple of the symbol freq.

% sampleTime = 20.5e-9;
sampleFreq = 100e6;
sampleTime = 1/sampleFreq;
desChirpTime = 101.1e-6;
desSymFreq = 1e6;

rollOff = 0.5;
span = 2; % This must be an even number

redPhase = pi/12; % 20 degrees

f0 = 4e6;
f1 = 12e6;
chirpPhase = 0;
Ac = 1;

[sentSignal, timeVec, RRCtransSignal, complexChirp, freqSymbol] = BRPSKmodChirpSignalGeneration(sampleTime, desChirpTime, desSymFreq, rollOff, span, redPhase, f0, f1, chirpPhase, Ac);

% figure;
% plot(timeVec, sentSignal);
% grid on;

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
figure;
plot(alphaFreq, db(mean(abs((SCF)), 2)));
grid on;