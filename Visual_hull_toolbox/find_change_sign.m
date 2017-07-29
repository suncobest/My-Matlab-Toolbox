function p0 = find_change_sign( xx )
%FIND_CHANGE_SIGN find where sign of vector XX change from minus to plus or vice versa. 
%   xx: input vector;
%   p0: the exact position of zero in unit of float index.
%   See also find_minus2plus, find_plus2minus, find_localmin, find_localmax,  find_localextreme.

%   ind: index of xx where the value is minus, but next element is 0 or positive,
%          and where the value is positive, but next element is 0 or negative;
assert(isvector(xx), 'Input must be a vector!');
ind = find((xx(2:end)>=0 & xx(1:end-1)<0) | (xx(2:end)<=0 & xx(1:end-1)>0));
if isempty(ind),
    p0 = [];
else
    k = xx(ind+1)-xx(ind);  % k~=0;
    p0 = ind-xx(ind)./k;      % x0 = x1-y1/k;
end;

return;


%% test
a = randn(1,10);
a = spline_interp3(1:10,a,0:0.1:10);
p1 = find_minus2plus(a);
p2 = find_plus2minus(a);
p = find_change_sign(a);