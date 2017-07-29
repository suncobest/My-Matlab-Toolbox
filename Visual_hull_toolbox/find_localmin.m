function [pmin, xmin] = find_localmin( xx )
%FIND_LOCALMIN find the stagnation points where vector XX achieve its loacal minimum.
%   xx: input vector;
%   pmin: the exact position of local minimum in unit of float index.
%   xmin: local minimum value of xx where derivative change sign;
%   See also find_localmax,  find_localextreme, find_minus2plus,
%   find_plus2minus, find_change_sign.

assert(isvector(xx), 'Input must be a vector!');
if length(xx)<3,
    pmin = [];
    xmin = [];
    return;
end;
% centered difference
dxx = xx(2:end)-xx(1:end-1);   % length=n-1, position of dxx (1:n-1)+0.5
ind = find(dxx(2:end)>=0 & dxx(1:end-1)<0);
if isempty(ind),
    pmin = [];
    xmin = [];
else
    k = dxx(ind+1)-dxx(ind);  % k>0;
    pmin = ind-dxx(ind)./k+0.5;      % x0 = x1-y1/k;
    ind = round(pmin);
    xmin = xx(ind);
end;

return;


%% test
a = rand(1,10)*10;
a = spline_interp3(1:10, a, (1:0.1:10)');
[p1, a1] = find_localmin(a);
[p2, a2] = find_localmax(a);