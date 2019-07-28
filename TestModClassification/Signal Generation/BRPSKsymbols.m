function modSymbols = BRPSKsymbols(binArray, redPhase)
% 
% modSymbols = modBRPSKsymbols(binArray, redPhase)
% 
% Computes the modulated BRPSK symbols from the binary input symbols
% 
% INPUTS:
% binArray      - binary array data
% redPhase      - reduced phase
% 
% OUTPUTS:
% modSymbols    - modulated BRPSK symbols
% 

%% Introductino of the reduced phase offset
introducedPhase = binArray*redPhase;

%% BRPSK modulation
modSymbols = exp(1i*introducedPhase);
end