function [l, u, id, detm, x] = bandlu(A, m1, m2, b)
% BANDLU - LU decompose the compact form of band-diagonal matrix A, and
% compute the determinant of A. Solve the equation band2mat(A,m1,m2)*X=b,
% when b exsit and is not empty.
%
% Expression: [L, U, ID, DETM, X] = BANDLU(IN, M1, M2, B).
% Algorithm: Numerical.Recipes.C++.3rd.Edition, P59.
%
% Given an n*m compactly stored band-diagonal matrix A with M1 subdiagonal
% rows and M2 superdiagonal rows, an LU decomposition of a rowwise
% permutation of A is constructed. The upper and lower triangular matrices
% are stored in U and L respectively.  The vector INDEX records the row
% permutation effected by the partial pivoting; The diagonal elements of U
% (whose product, times d=+-1, gives the determinant DETM) are in the first
% column of U. (d is +-1 depending on whether the number of row interchanges
% was even or odd respectively.)
%
% ID is more important than INDEX. ID records the process of permutation
% and can lead to the result INDEX.
%
% When B exist and is not empty, X is the solution of  band2mat(A,m1,m2)*X=B.
% See also ludcomp, lu, bandmul, band_transpose.

% by zpf, form BIT, 2015-7-6

assert(ismatrix(A) && ~isempty(A), 'The 1st input argument must be a non-empty matrix!');

[n,m] = size(A);       % n rows, m cols
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

assert(min(m1,m2)>=0 && max(m1,m2)<=n-1, ['The 2nd and 3rd input arguments ' ...
                    'must be non-negative and fit the matrix size!']);

if nargin == 4,
    assert(ismatrix(b) && ~isempty(b), 'The 4th input argument must be a non-empty matrix!');
    [nx, ~] = size(b);
    assert(nx==n,'Number of rows of the 1st and 4th input arguments must be equal!');
end;

TINY = 1.0e-40;
id  = zeros(1,n);
%   index = 1:n;

% Rearrange A: to drive out zeros in the first columnn
d = A(1:m1, :);
A(1:m1, :) = 0;
for i=1:m1,
    A(i,1:m2+i) = d(i,m1+2-i:m);
end;

d = 1;
l = zeros(n,m1);
for k=1:n,
    tmpend = min(m1+k,n);
    [pivot, i] = max(abs(A(k:tmpend,1)));    % Find the pivot element.
    i = i+k-1;
    id(k) = i;
    if k<i,
        d = -d;
        %   index([k,i]) = index([i,k]);
        A([k,i],:) = A([i,k],:);        % only u is interchanged, l must be transformed after
    end;
    % Matrix is algorithmically singular, but proceed anyway with TINY pivot.
    if pivot==0,
        warning('The input matrix is singular!');
        A(k,1)=TINY;
    end;

    % Do the elimination.
    dum = A(k+1:tmpend,1)/A(k,1);
    l(k, 1:tmpend-k) = dum;
    A(k+1:tmpend, 1:m-1) = A(k+1:tmpend, 2:m) - dum*A(k, 2:m);
    A(k+1:tmpend, m) = 0;
end;

tmpend = min(m,n);
u = A(:,1:tmpend);      % the compact form of upper tridiagonal matrix should not exceed n column
detm = d*prod(u(:,1));

if nargin<4,
    x = [];
    return;
end;

x = b;
for k=1:n,
    i = id(k);
    if i~=k,
        % interchange x like transforming l to lower tridiagonal matrix
        x([i, k],:) = x([k, i],:);
    end;
    tmpend = min(m1+k, n);
    A = l(k, 1:tmpend-k)';
    x(k+1:tmpend, :) = x(k+1:tmpend, :) - A*x(k, :);   % forward substitution
end;

x(n,:) = x(n,:)/u(n,1);
for k=n-1:-1:1,
    tmpend = min(m, n+1-k);
    x(k,:) = (x(k,:)-u(k, 2:tmpend)*x(k+1:tmpend+k-1,:))/u(k,1);    %  backsubstitution
end;


return;



%% transform l to lower tridiagonal matrix 'lower'
lower = zeros(n);
index = 1:n;
for k = 1:n,
    i = id(k);
    if i~=k,
        index([i, k]) = index([k, i]);
        lower([i, k],1:k-1) = lower([k, i],1:k-1);
    end;
    tmpend = min(k+m1, n);
    lower(k+1:tmpend, k) = l(k, 1:tmpend-k);
end;
lower = lower + eye(n);


%% Test
n = 6;
m1 = 5;
m2 = 2;
am = 10*randn(n);
ab = mat2band(am, m1, m2);
am = band2mat(ab, m1);
b = 10*randn(n,m1);
[l, u, id, detm, x] = bandlu(ab, m1, [], b);
[ll, uu, p] = lu(am,'vector');

lower = zeros(n);
index = 1:n;
for k = 1:n,
    i = id(k);
    if i~=k,
        index([i, k]) = index([k, i]);
        lower([i, k],1:k-1) = lower([k, i],1:k-1);
    end;
    tmpend = min(k+m1, n);
    lower(k+1:tmpend, k) = l(k, 1:tmpend-k);
end;
lower = lower + eye(n);
err = am(index,:)-lower*band2mat(u,0)
am*x-b

%%      solve x*am =b, so am.'*x.'=b.'
n=6;
m1 = 5;
m2 = 2;
am = 10*randn(n);
ab = mat2band(am, m1, m2);
am = band2mat(ab, m1);
b = 10*randn(m1,n);
[abt, mm1, mm2] = band_transpose(ab,m1);
[l, u, id, detm, x] = bandlu(abt, mm1, [], b');
x'*am-b