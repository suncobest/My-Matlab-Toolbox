function [pex, xex, idmin] = find_localextreme( xx )
%FIND_LOCALEXTREME find the stagnation points where vector XX achieve its loacal extreme value.
%   xx: input vector;
%   pex: the exact position of local extreme in unit of float index.
%   xex: local extreme value of xx where derivative change sign;
%   idmin: logical index of minimum points in xex;
%        then the logical index of maximum points in xex is ~idmin;  
%   See also find_localmin, find_localmax, find_minus2plus,
%   find_plus2minus, find_change_sign.

assert(isvector(xx), 'Input must be a vector!');
if length(xx)<3,
    pex = [];
    xex = [];
    idmin = [];
    return;
end;
% centered difference
dxx = xx(2:end)-xx(1:end-1);   % length=n-1, position of dxx (1:n-1)+0.5
ind = find((dxx(2:end)>=0 & dxx(1:end-1)<0) | (dxx(2:end)<=0 & dxx(1:end-1)>0));
if isempty(ind),
    pex = [];
    xex = [];
    idmin = [];
else
    k = dxx(ind+1)-dxx(ind);  % k~=0;
    pex = ind-dxx(ind)./k+0.5;      % x0 = x1-y1/k;
    ind = round(pex);
    xex = xx(ind);
    idmin = dxx(ind)-dxx(ind-1)>0;  % 2nd derivative > 0 to reach local minimum
end;

return;



%% test
a = rand(1,10)*10;
a = spline_interp3(1:10, a, 1:0.1:10);
[p0, a0, id] = find_localextreme(a);
