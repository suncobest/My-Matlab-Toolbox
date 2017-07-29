function p0 = find_plus2minus( xx )
%FIND_PLUS2MINUS find where sign of vector XX change from plus to minus. 
%   xx: input vector;
%   p0: the exact position of zero in unit of float index.
%   See also find_minus2plus, find_change_sign, find_localmin,
%   find_localmax,  find_localextreme.

%   ind: index of xx where the value is positive, but next element is 0 or negative;
assert(isvector(xx), 'Input must be a vector!');
ind = find(xx(2:end)<=0 & xx(1:end-1)>0);
if isempty(ind),
    p0 = [];
else
    k = xx(ind+1)-xx(ind);  % k<0;
    p0 = ind-xx(ind)./k;      % x0 = x1-y1/k;
end;

return;

