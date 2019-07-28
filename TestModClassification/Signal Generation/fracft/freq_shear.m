function y = freq_shear(x, c)
% freq_shear -- chirp multiplication
%
%  Usage
%    y = freq_shear(x, c)
%
%  Inputs
%    x     signal vector, must have an odd length (to have a center point)
%    c     chirp rate
%
%  Outputs
%    y     the frequency sheared signal
%
% If x has length N and bandwidth B, then y will have length N and 
% bandwidth c*N*B/pi.  One must be careful to avoid aliasing.

% Copyright (C) -- see DiscreteTFDs/Copyright

error(nargchk(2, 2, nargin));

x = x(:);
N = length(x);
if (rem(N,2)==0)
  error('signal length must be odd')
end
M = (N-1)/2;
n = (-M:M)';

y = x .* exp(1j*c/2*n.^2);
