function [out, mm1, mm2] = band_transpose(in, m1, m2)
% BAND_TRANSPOSE - transpose a compact format band-diagonal matrix.
%
%  Expression: [OUT, MM1, MM2] = BAND_TRANSPOSE(IN, M1, M2).
%
%  IN is a compact format band-diagonal matrix, with N rows and M1+M2+1
% columns. M1 denotes the band width of subdiagonal elements, and M2 the
% width of superdiagonal ones. OUT is the transposed compact format
% band-diagonal matrix.
%
% Tilt the compact format 45 degree clockwise in the y direction, then
% fliplr the matrix and swap the value of M1 and M2 to get results. Some
% subdiagonal and superdiagonal elements that are outside the square matrix
% will be seem as zeros. See also band2mat, mat2band, bandmul.

% by zpf, form BIT, 2015-7-6

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

out = zeros(n, m);
for i =1:m,
    k = i-m1-1;
    tmpstart = max(1,1-k);
    tmpend = min(n, n-k);
    out((tmpstart:tmpend)+k, i) = in(tmpstart:tmpend, i);
end;
out = out(:, end:-1:1);
mm1 = m2;
mm2 = m1;

return;



%% Test
a = randi(10,6,6)
a = mat2band(a, 2)
[b, m1, m2] = band_transpose(a, 2)

%%
a = randi(10,6,11)
b = band_transpose(a)
[a, m1, m2] = band_transpose(b)