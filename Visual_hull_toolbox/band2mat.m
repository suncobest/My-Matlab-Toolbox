function [out] = band2mat(in, m1, m2)
% BAND2MAT - transform a compact format band-diagonal matrix to a full
% square one.
%  Expression: OUT = BAND2MAT(IN, M1, M2).
%
%  IN is a compact format band-diagonal matrix, with N rows and M1+M2+1
% columns. M1 denotes the band width of subdiagonal elements, and M2 the
% width of superdiagonal ones. OUT is the output full squre band-diagonal
% matrix.
%
% Tilt the compact format 45 degree counter-clockwise in the x direction,
% then the full one is there. Some subdiagonal and superdiagonal elements
% will be squeezed out.
% See also mat2band.

% by zpf, form BIT, 2015-7-2

assert(ismatrix(in),'The 1st input argument must be a matrix!');
if isempty(in),
    out = [];
    return;
end;

[n,m] = size(in);       % n rows, m cols
assert(m<=2*n-1,'The band-diagonal width should not exceed the matrix size!');
if nargin<3 || isempty(m2),
    if nargin<2 || isempty(m1),
        m1 = n-1;
        m2 = m-m1-1;
    else
        m1 = floor(m1);
        m2 = m-m1-1;
    end;
else
    if nargin<2 || isempty(m1),
        m2 = floor(m2);
        m1 = m-m2-1;
    else
        m1 = floor(m1);
        m2 = floor(m2);
        assert(m1+m2+1==m,'The band-diagonal width do not match the input compact form!');
    end;
end;

assert(min(m1,m2)>=0 && max(m1,m2)<=n-1,['The 2nd and 3rd input arguments ' ...
                    'must be non-negative and fit the matrix size!']);

out = zeros(n);
% algorithm: Numerical.Recipes.C++.3rd.Edition, P59.
for i=1:n,
    % k is the distance from the diagonal column (m1+1) to i th column in the
    % compact form. Getting a minus triangle and a positive one straddle
    % the diagonal column, like the stress distribution of bending moment.
    k = i-m1-1;
    tmpstart = max(1,1-k);
    tmpend = min(m, n-k);
    out(i,(tmpstart:tmpend)+k) = in(i,tmpstart:tmpend);
end;

return;


if 0,           % equivalent algorithm as above
    for i=1:m1,
        a = min(n,m2+i);
        b = m1+2-i;
        c = a-1+b;
        out(i,1:a) = in(i,b:c);
    end;
    for i=m1+1:n-m2,   % work when m<=n (m<=m2+i<=n, so min(n,m2+i)=m2+i)
        out(i,i-m1:m2+i) = in(i,1:m);
    end;
    a = max(m1+1,n-m2+1);
    for i = a:n,
        out(i,i-m1:n) = in(i,1:n+m1+1-i);
    end;
end;




%% Test
a = randi(10,6,6)
b = band2mat(a,[],2)

%%
a = randi(10,6,11)
b = band2mat(a)