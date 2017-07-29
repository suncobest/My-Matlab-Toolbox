function [out] = mat2band(in, m1, m2)
% MAT2BAND - transform a full square matrix to a compact format
% band-diagonal matrix.  Expression: OUT = MAT2BAND(IN, M1, M2).
%
%  IN is a full squre matrix, assumed to be size of N^2. M1 denotes the band
% width of subdiagonal elements, and M2 the width of superdiagonal ones. OUT
% is the output compact format band-diagonal matrix of N rows and M1+M2+1
% columns.
%
% Tilt the full square matrix 45 degree clockwise in the x direction, then
% you got the compact format band-diagonal matrix. Subdiagonal and
% superdiagonal elements that are outside the square matrix will be seem as
% zeros.  See also band2mat.

% by zpf, form BIT, 2015-7-2

assert(ismatrix(in),'The 1st input argument must be a matrix!');
if isempty(in),
    out = [];
    return;
end;

[n,m] = size(in);
assert(n==m,'The input matrix must be square matrix!');
if nargin<3 || isempty(m2),
    if nargin<2 || isempty(m1),
        m1 = n-1;
        m2 = n-1;
    else
        m1 = floor(m1);
        m2 = n-1;
    end;
else
    if nargin<2 || isempty(m1),
        m1 = n-1;
        m2 = floor(m2);
    else
        m1 = floor(m1);
        m2 = floor(m2);
    end;
end;

assert(min(m1,m2)>=0 && max(m1,m2)<=n-1,['The 2nd and 3rd input arguments ' ...
                    'must be non-negative and fit the matrix size!']);

m = m1+m2+1;
out = zeros(n,m);

% algorithm: Numerical.Recipes.C++.3rd.Edition, P59.
for i=1:n,
    % k is the distance from the diagonal column (m1+1) to i th column in the
    % compact form. Getting a minus triangle and a positive one straddle
    % the diagonal column, like the stress distribution of bending moment.
    k = i-m1-1;         
    tmpstart = max(1,1-k);
    tmpend = min(m, n-k);
    out(i,tmpstart:tmpend) = in(i,(tmpstart:tmpend)+k);
end;

return;


if 0,           % equivalent algorithm as above
    for i=1:m1,
        a = min(n,m2+i);
        b = m1+2-i;
        c = a-1+b;
        out(i,b:c) = in(i,1:a);
    end;
    for i=m1+1:n-m2,   % work when m<=n (m<=m2+i<=n, so min(n,m2+i)=m2+i)
        out(i,1:m) = in(i,i-m1:m2+i);
    end;
    a = max(m1+1,n-m2+1);
    for i = a:n,
        out(i,1:n+m1+1-i) = in(i,i-m1:n);
    end;
end;




%% Test
a = randi(10,6,6)
b = mat2band(a,[],2)

%%
a = randi(10,5,5)
b = band2mat(mat2band(a))


