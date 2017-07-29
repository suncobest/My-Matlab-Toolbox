% below result is only accurate at decimal 0.1, but not the tolerate, why?
% is the conjugate gradient algorithm only for square matrix A's system?
% about the lambda1 + lambda2 = 1, how embedded this in algorithm?
% below is Matlab Code

w = CGResult(1);
lambda1 = CGResult(2);
lambda2 = CGResult(3);
(w+13/35*lambda1+47/140*lambda2) >= 2/5
(w+57/140*lambda1+61/140*lambda2) >= 2/5
(w+31/140*lambda1+8/35*lambda2) >= 1/5


A = [1, 13/35, 47/140;
1, 57/140, 61/140;
1, 31/140, 8/35;];
b = [2/5;2/5;1/5;];
tol = 0.0000000000001;
Max = 10000;
x2 = [0.5;0.5;0.5];
n2 = 3;
MartinCG(x2, A, b, Max, tol, n2);

% Why only for Square Matrix, How about Non Square Matrix A

function CGResult = MartinCG(x2, A, b, Max, tol, n3)
r = zeros(n3);
x = zeros(n3,Max-1);
x(:,1) = x2;
r = zeros(n3,Max-1);
v = zeros(n3,Max-1);
t = zeros(n3,Max-1);
dum = zeros(n,n);
r(:,1) = b - A*x(:,1);
v(:,1) = r(:,1);
for k=1:1:Max-3
    k
    testnorm = norm(v(:,k))
    if (norm(v(:,k)) == 0)
        break;
    end
    t(:,k) = (r(:,k).*r(:,k))./(v(:,k).*(A*v(:,k)));
    x(:,k+1) = x(:,k) + (t(:,k).*v(:,k));
    r(:,k+1) = r(:,k) - t(:,k).*(A*v(:,k));
    testnorm2 = norm(r(:,k+1))
    if (norm(r(:,k+1)) < tol)
        break;
    end
    v(:,k+1) = r(:,k+1) + (r(:,k+1).*r(:,k+1))./(r(:,k).*r(:,k)).*v(:,k);
end
CGResult = x(:,k+1);
end
