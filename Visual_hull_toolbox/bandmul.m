function y = bandmul(A, m1, m2, x)
% BANDMUL - Band-diagonal matrix A left multiply another matrix X.
%   Expression: Y = BANDMUL(A,M1,M2,X) is equivalent to Y =
%   BAND2MAT(A,M1,M2)*X, but solved more efficiently. Only the nonzero
%   elements of full square format of matrix A is chosen for caculation.
%
%   A is a compact format band-diagonal matrix, with N rows and M1+M2+1
%   columns. M1 denotes the band width of subdiagonal elements, and M2 the
%   width of superdiagonal ones.
%   See also bandmul2, bandlu, band_transpose.

%   by zpf, form BIT, 2015-7-3.

assert(ismatrix(A) && ismatrix(x),'The 1st input argument must be a matrix!');
if isempty(A) || isempty(x),
    y = [];
    return;
end;

[n,m] = size(A);       % n rows, m cols
assert(m<=2*n-1,'The band-diagonal width should not exceed the matrix size!');

[nx, mx] = size(x);
assert(nx==n,'Number of rows of the 1st and 4th arguments must be equal!');

m1 = floor(m1);
m2 = floor(m2);

if isempty(m2),
    if isempty(m1),
        m1 = 0;
        m2 = m-1;
    else
        m2 = m-m1-1;
    end;
else
    if isempty(m1),
        m1 = m-m2-1;
    else
        assert(m1+m2+1==m,'The band-diagonal width do not match the input compact form!');
    end;
end;

assert(min(m1,m2)>=0 && max(m1,m2)<=n-1,['The 2nd and 3rd input arguments ' ...
    'must be non-negative and fit the matrix size!']);

y = zeros(n,mx);

% algorithm: Numerical.Recipes.C++.3rd.Edition, P59.
for i=1:n,
    % k is the distance from the diagonal column (m1+1) to i th column in the
    % compact form. Getting a minus triangle and a positive one straddle
    % the diagonal column, like the stress distribution of bending moment.
    k = i-m1-1;
    tmpstart = max(1,1-k);
    tmpend = min(m, n-k);
    y(i,1:mx) = A(i,tmpstart:tmpend)*x((tmpstart:tmpend)+k,1:mx);
end;

return;


%% Test
a = magic(6)
a(:,1:2)=[]
c = randi(10,6,3)
b = bandmul(a,[],[],c)

%%
a = randi(10,6,10)
b = bandmul(a,4,[],c)
