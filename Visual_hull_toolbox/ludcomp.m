function [l, u, index, detm, x] = ludcomp(A, b)
% LUDCOMP - Given a N*N matrix A, this routine return the LU decomposition
% of a rowwise permutation of itself. Solve the equation A*X=B, when B exsit
% and is not empty.
%
%   Algorithm: Numerical.Recipes.C++.3rd.Edition, P52.
%
%   LU = LUDCOMP(A). A is input. On output, LU contains the strict lower
%   triangle of L embedded in the same matrix as the upper triangle of
%   U. The permutation information is lost.
%
%   [L, U] = LUDCOMP(A) stores an upper triangular matrix in U and a
%   "psychologically lower triangular matrix" (i.e. a product of lower
%   triangular and permutation matrices) in L, so that A = L*U.
%
%   [L, U, INDEX] = LUDCOMP(A).  INDEX is an output vector that records the
%   row permutation effected by the partial pivoting; so that A(INDEX,:) =
%   L*U. This routine is used in combination with solve to solve linear
%   equations or invert a matrix.
%
%    [L, U, INDEX, DETM] = LUDCOMP(A).  DETM is the determinant of matrix A.
%
%    [L, U, INDEX, DETM, X] = LUDCOMP(A, B).  B is input right hand, X is
%    the solution of A*X=B.
%   See also lu, bandlu.

% by zpf, form BIT, 2015-7-7

assert(ismatrix(A) && ~isempty(A), 'The 1st input argument must be a non-empty matrix!');

[n,m] = size(A);
assert(n==m,'The input matrix must be square matrix!');

if nargin == 2,
    assert(ismatrix(b) && ~isempty(b),['The 2nd input argument must be a non-empty matrix!']);
    [nx, ~] = size(b);
    assert(nx==n,'Number of rows of the two input arguments must be equal!');
end;

vv = max(abs(A),[],2);
if all(vv)==0,
    error('The input matrix is singular, some rows are all zeros!');
end;

% Implicit pivoting: vv stores the implicit scaling of each row.
A = A./vv(:,ones(1,n));
TINY = 1.0e-40;
index = 1:n;
%   id  = zeros(1,n);

d = 1;
for k=1:n,
    [pivot, i] = max(abs(A(k:n,k)));    % Find the pivot element.
    %   id(k) = i+k-1;
    if i>1,
        %   d is 1 or -1 depending on whether the number of row
        %   interchanges is even or odd respectively.
        d = -d;
        i = i+k-1;
        index([k,i]) = index([i,k]);
        A([k,i],:) = A([i,k],:);        % interchange l and u at the same time
        vv([k,i]) = vv([i,k]);
    end;
    % Matrix is algorithmically singular, but proceed anyway with TINY pivot.
    if pivot==0,
        warning('The input matrix is singular!');
        A(k,k)=TINY;
    end;

    % Do the elimination.
    dum = A(k+1:n,k)/A(k,k);
    A(k+1:n, k) = dum;
    A(k+1:n, k+1:n) = A(k+1:n, k+1:n) - dum*A(k, k+1:n);
end;

A = A.*vv(:,ones(1,n));         %   inv(S)*A = L*U, A = S*L*inv(S)*S*U.
for i=1:n-1,
    A(i+1:n, i)=A(i+1:n, i)/vv(i);
end;

l = A;
if nargout <2,
    return;
end;

for i=1:n,
    l(i,i) = 1;
    l(i, i+1:n)=0;
end;

u = A;
detm = d;
for i=1:n,
    u(i,1:i-1)=0;
    detm = detm*u(i,i);
end;

if nargout < 3,
    [~, A] = sort(index);
    l = l(A,:);
    return;
end;

if nargin < 2,
    x = [];
    return;
end;

x = b;
x = x(index,:);
for i=2:n,
    x(i,:) = x(i,:)-A(i,1:i-1)*x(1:i-1,:);    % forward substitution
end;

x(n,:) = x(n,:)/A(n,n);
for i=n-1:-1:1,
    x(i,:) = (x(i,:)-A(i, i+1:n)*x(i+1:n,:))/A(i,i);    %  backsubstitution
end;

return;



%% Test
n =6;
a = 10*rand(n);
[lower, upper, p] = lu(a,'vector');
[l, u, index, detm] = ludcomp(a);
err1= a(index,:) - l*u
err2 = a(p,:) -lower*upper

%%
a = 10*randn(8);
b = 10*randn(8,3);
x = a\b;
[~,~,~,d,xx] = ludcomp(a,b);
x-xx
