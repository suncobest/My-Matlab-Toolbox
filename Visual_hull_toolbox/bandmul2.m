function [y1,y2,t1,t2] = bandmul2(A, m1, m2, x)
% BANDMUL2 - Band-diagonal matrix A left multiply another matrix X.
%   Expression: Y = BANDMUL2(A,M1,M2,X) is equivalent to Y =
%   BAND2MAT(A,M1,M2)*X, but solved more efficiently. Only the nonzero
%   elements of full square format of matrix A is chosen for caculation.
%
%   A is a compact format band-diagonal matrix, with N rows and M1+M2+1
%   columns. M1 denotes the band width of subdiagonal elements, and M2 the
%   width of superdiagonal ones.
%   Compare direct matrix multiplication and loops. 
%   See also banddiv.

%   by zpf, form BIT, 2015-7-3.

assert(ismatrix(A) && ismatrix(x),'The 1st input argument must be a matrix!');
if isempty(A) || isempty(x),
    y1 = [];
    return;
end;

[n,m] = size(A);       % n rows, m cols
assert(m<=2*n-1,'The band-diagonal width should not exceed the matrix size!');

[nx, mx] = size(x);
assert(nx==n,'The rows of the 1st and 4th arguments must be equal!');

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

y1 = zeros(n,mx);
tic;
for i=1:n,
    k = i-m1-1;
    tmpstart = max(1,1-k);
    tmpend = min(m, n-k);
    y1(i,1:mx) = A(i,tmpstart:tmpend)*x((tmpstart:tmpend)+k,1:mx);
end;
t1 = toc;

% return;

y2 = zeros(n,mx);
tic;
% algorithm: Numerical.Recipes.C++.3rd.Edition, P59.
% equivalent algorithm as above
for i=1:n,
    k = i-m1-1;
    tmpstart = max(1,1-k);
    tmplend = min(m, n-k);
    for nl =1:mx,
        for j = tmpstart:tmplend,
            y2(i,nl) = y2(i,nl)+A(i,j)*x(j+k,nl);
        end;
    end;
end;
t2 = toc;

return;


%% Test
% loops are faster when matrix size is small, like less than 200.
% when the scale goes up, direct matrix manipulation is getting more efficient. 
n1 = 265;          % 20, 1000
a = randn(n1);
c = randi(10,n1,4);
[b1,b2,t1,t2] = bandmul2(a,4,[],c);
norm(b1(:)-b2(:))
t1-t2