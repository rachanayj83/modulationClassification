function y = time_shear(x, c)
% time_shear -- chirp convolution
%
%  Usage
%    y = time_shear(x, c)
%
%  Inputs
%    x     signal vector, must have an odd length (to have a center point)
%    c     chirp rate
%
%  Outputs
%    y     the time sheared signal
%
% If x has bandwidth B and duration d, then y will have bandwidth B and 
% duration c*N*d/pi.  One must be careful to avoid aliasing in time.

% Copyright (C) -- see DiscreteTFDs/Copyright

error(nargchk(2, 2, nargin));

x = x(:);
N = length(x);
if (rem(N,2)==0)
  error('signal length must be odd')
end
M = (N-1)/2;

interp = ceil(2*abs(c)/(2*pi/N));
xx = sinc_interp(x,interp)/interp;
n = (-2*M:1/interp:2*M)';
z = exp(1j*c/2*n.^2);
y = lconv(xx, z);
center = (length(y)+1)/2;
y = y(center-interp*M:interp:center+interp*M);
y = y*sqrt(c/2/pi);

% bad way -- inaccurate answers for large c
%interp = 2;
%xx = sinc_interp(x,interp);
%n = [-4*M*(2*pi/N)/c:4*M*(2*pi/N)/c]'/2;
%z = exp(j*c/2*n.^2);
%y = lconv(xx, z);
%center = (length(y)+1)/2;
%y = y(center-interp*M:interp:center+interp*M);

%%% another bad way -- does circular convolutions
%center = (N+1)/2;
%x = [x(center:end) ; x(1:center-1)];
%X = fft(x);
%X = [X(center+1:end) ; X(1:center)];
%Y = time_shear(X,-c);
%Y = [Y(center:end) ; Y(1:center-1)];
%y = ifft(Y);
%y = [y(center+1:end) ; y(1:center)];

