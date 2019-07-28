function y = fracft(x,a)
% fracft -- compute the fractional Fourier transform
%
%  Usage
%    y = fracft(x,a)
%
%  Inputs
%    x     signal vector, must have an odd length (to have a center point)
%    a     fraction. a=1 corresponds to the Fourier transform and a=4
%          is an identity transform.
%
%  Outputs
%    y     the fractional Fourier transform
%
% Example:
%   x = chirplets(63,[1 48 pi/2 0 sqrt(63/4/pi)]);
%   for i=0:0.25:4,
%     y = fracft(x,i); wigner1(y); axis square; pause 
%   end
%
% Algorithm taken from H. M. Ozaktas et al., Digital Computation of the 
% Fractional Fourier Transform, IEEE Trans. Signal Processing, September,
% 1996.

% Copyright (C) -- see DiscreteTFDs/Copyright

error(nargchk(2, 2, nargin));

x = x(:);
N = length(x);
M = (N-1)/2;

if (rem(N,2)==0)
  error('signal length must be odd')
end

% do special cases
a = mod(a,4);
if (a==0)
  y = x;
  return
elseif (a==1)
  y = fft([x(M+1:N); x(1:M)])/sqrt(N);
  y = [y(M+2:N); y(1:M+1)];
  return
elseif (a==2)
  y = flipud(x);
  return
elseif (a==3)
  y = ifft([x(M+1:N); x(1:M)])*sqrt(N);
  y = [y(M+2:N); y(1:M+1)];
  return
end

% get a in the right range
if (a>2)
  x = flipud(x);
  a = a-2;
end
if (a>1.5)
  x = fft([x(M+1:N); x(1:M)])/sqrt(N);
  x = [x(M+2:N); x(1:M+1)];
  a = a - 1;
end
if (a<0.5)
  x = ifft([x(M+1:N); x(1:M)])*sqrt(N);
  x = [x(M+2:N); x(1:M+1)];
  a = a + 1;
end

phi = a*pi/2;
alpha = cot(phi);
beta = csc(phi);

P = 3;  % oversampling factor
x = sinc_interp(x,P);
x = [zeros(2*M,1) ; x ; zeros(2*M,1)];

x2 = freq_shear(x, 2*pi/N*(alpha-beta)/P^2);
x3 = time_shear(x2, 2*pi/N*beta/P^2);
y = freq_shear(x3, 2*pi/N*(alpha-beta)/P^2);

y = y(2*M+1:end-2*M);
y = y(1:P:end);
y = exp(-1j*(pi*sign(sin(phi))/4-phi/2))*y;
