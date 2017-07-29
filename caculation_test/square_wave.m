function [f,c] = square_wave(A,L,N)
%SQUARE_WAVE use Fourier series to get an approximation.
%   The function is odd, so only the SIN coefficient left. A is the
%   amplitude. 2L is the period. N is the highest order of coefficient. The
%   square wave f(x) is -A between [-L,0) and A between [0, L).

c=linspace(-2*L,2*L,500);   % x sample
r=1:2:N;    % only odd Sin coefficient left.
[C,R]=meshgrid(c,r);
x=R.*C;
if isvector(x)
    f=4*A/pi*((1./R).*sin(pi/L*x));
else
    f=sum(4*A/pi*((1./R).*sin(pi/L*x)));
end

end

