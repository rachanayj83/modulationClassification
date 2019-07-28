function y = lconv(x, h)
% lconv -- perform a linear convolution with ffts
%
%  Usage
%    y = lconv(x, h)
%
%  Inputs
%    x, h   input vectors
%
%  Outputs
%    y      the linear convolution of x and h

% Copyright (C) -- see DiscreteTFDs/Copyright

error(nargchk(2, 2, nargin));

x = x(:);
N = length(x);
h = h(:);
M = length(h);
P = 2^nextpow2(N+M-1);

y = ifft( fft(x,P) .* fft(h,P));
y = y(1:N+M-1);

if (isreal(x) & isreal(h))
  y = real(y);
end
