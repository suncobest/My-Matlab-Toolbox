function [mxx, stdxx] = period_average(xx, nc)
%PERIOD_AVERAGE compute average of XX over all periods.
%   XX: periodically changed value. (m*n), m is dimension of the value.
%   NC: denotes N columns in 1 period.
%   MXX: mean of each dimension of XX over all cycles.
%   STDXX: standard deviation of each dimension of XX over all cycles.

assert(ismatrix(xx),'The 1st argument must be matrix!');
[m,n] = size(xx);
np = n/nc;
assert(np==round(np),'Number of columns of the 1st argument must be divided exactly by the 2nd input!');
xx = reshape(xx,[m,nc,np]);
mxx = mean(xx,3);
xx = (xx-mxx(:,:,ones(1,np))).^2;
stdxx = sqrt(sum(xx,3)/(np-1));

return;


%% test
x = rand(1,100);
[mx, dx] = period_average(x, 10);
x = reshape(x,10,[])';
mx1 = mean(x);
dx1 = std(x);
norm(mx-mx1)
norm(dx-dx1)