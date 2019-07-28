function out = sinc_interp(x, a)
% sinc_interp -- sinc interpolate a signal
%
%  Usage
%    y = sinc_interp(x, a)
%
%  Inputs
%    x      signal vector
%    a      interpolation factor (optional, default is 2)
%
%  Outputs
%    y     interpolated vector.  If N=length(x) then
%          length(y) = a*N-a+1 and y(1:a:end) = x.
%          No extrapolation is done.
%
% Surprisingly this does not seem to be included with matlab.

% Copyright (C) -- see DiscreteTFDs/Copyright

error(nargchk(1, 2, nargin));

if (nargin < 2)
  a = 2;
end

x = x(:);
N = length(x);
M = a*N-a+1;

% y has length: a*N-a+1
y = zeros(M,1);
y(1:a:M) = x;

% h has length: 2*(a*N-a-1)+1
h = sinc([-(N-1-1/a):1/a:(N-1-1/a)]');

% out has length 3*(a*N-a)-1
out = lconv(y, h);

% what we want has length: a*N-a+1
out = out(a*N-a:end-a*N+a+1);
