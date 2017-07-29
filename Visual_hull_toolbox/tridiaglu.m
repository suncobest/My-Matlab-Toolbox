function [l, u, detm, x] = tridiaglu(A, b, SW)
% TRIDIAGLU - lu decomposition of tridiagnal matrix.  LU decompose the
% compact form of tridiagonal matrix A (or transform a N*N full matrix to
% compact form --three vectors), and compute the determinant of A. Solve the
% equation band2mat(A,1,1)*X=B, when B exsit and is not empty (If A is a
% square matrix, then the equation is A*X =B).
% SW is the switch to turn on the function of decomposing full format matrix A
% if SW==1. If SW==0 or default, then the matrix A is treated as compact format.
%
% Expression: [L, U, DETM, X] = TRIDIAGLU(A,  B).
% Algorithm: Numerical.Recipes.C++.3rd.Edition, P56.
% See also bandlu, lu.

% by zpf, form BIT, 2015-7-1

assert(ismatrix(A) && ~isempty(A), 'The 1st input argument must be a non-empty matrix!');
[n, m]=size(A);
if nargin<3,
    SW = false;
else
    SW = ~~SW;
end;

if SW,
    assert(m==n && n>1,'The 1st matrix is taken as full format, so it must be square with dimension more than1!');
    va=[0; diag(A,-1)];
    vb=diag(A);
    vc=[diag(A,1); 0];
else
    assert(m==3 && n>=2,'The 1st matrix is taken as compact format matrix, so it must have 3 columns and at least two rows!');
    va = A(:,1);
    vb = A(:,2);
    vc = A(:,3);
end;

l = va;
u = vb;
bet = u(1);
detm = bet;
if nargin < 2,
    for i=2:n,
        if bet==0,
            error('Algorithm fails! Please recheck the matrix or try bandlu instead!');
        end;
        A = l(i)/bet;
        l(i) = A;
        bet=vb(i)-vc(i-1)*A;
        u(i) = bet;
        detm = detm*bet;
    end;
    x = [];
else
    assert(ismatrix(b) && ~isempty(b), 'The 2nd input argument must be a non-empty matrix!');
    [nx, ~] = size(b);
    assert(nx==n,'Number of rows of the two input arguments must be equal!');
    x = b;
    for i=2:n,
        if bet==0,
            error('Algorithm fails! Please recheck the matrix or try bandlu instead!');
        end;
        A = l(i)/bet;
        l(i) = A;
        x(i,:) = x(i,:) - A*x(i-1,:);
        bet=vb(i)-vc(i-1)*A;
        u(i) = bet;
        detm = detm*bet;
    end;
    if bet==0,
        error('Algorithm fails! The 1st matrix is singular!');
    end;
    x(n,:) = x(n,:)/bet;
    for i=n-1:-1:1,
        x(i,:) = (x(i,:)-vc(i)*x(i+1,:))/u(i);
    end;
end;

if nargout <=1,
    l=[l,u,vc];
    return;
end;

l=eye(n)+diag(l(2:n),-1);
u=diag(u)+diag(vc(1:n-1),1);

return;



%% test
n = randi(15);
a=diag(10*randn(n,1))+diag(10*randn(n-1,1),1)+diag(10*randn(n-1,1),-1);
b=randn(n, n-1);
[l,u,detm,x]=tridiaglu(a, b,1);
err1 = l*u-a
err2 = a*x-b

%%
n = randi(10);
ab=randn(n,3);
b=randn(n, n-1);
a=band2mat(ab,1);
[l,u,detm,x]=tridiaglu(ab,b);
err1 = l*u-a
err2 = a*x-b

