function binArray = binDataGeneration(numSymbols)
%
% binArray = binDataGeneration(numSymbols)
% 
% Generates a random binary, -1 and 1, data array. The array will have a 
% length of numSymbols.
%
% INPUT:
% numSymbols    - Number of symbols this will be the length of the array.
% 
% OUTPUT:
% binArray  - Array of binary data.
%

bitArray2 = unidrnd(2, [numSymbols, 1]); % We create an array of lenData values from a uniform distribution between 1 and 2
bitArray1 = bitArray2 - 1; % We substract 1 so that everything stays between 0 and 1
bitArray0 = bitArray2 -2; % We substract 2 so that everything stays between -1 and 0
binArray = bitArray1 + bitArray0; % We add the two arrays and we end up with binary data, an array with -1 and 1
end