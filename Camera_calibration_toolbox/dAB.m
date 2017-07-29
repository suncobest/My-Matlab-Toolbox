function [dABdA,dABdB] = dAB(A,B)

%      [dABdA,dABdB] = dAB(A,B);
%
%      returns : dABdA and dABdB

[p,n] = size(A); [n2,q] = size(B);

assert(n==n2,' A and B must have equal inner dimensions!');

nAB = p*q;
nA = p*n;
nB = q*n;

if issparse(A) ||  issparse(B) || nAB*nA>400 ,
  dABdA=sparse([],[],[],nAB,nA,p*nB);
else
  dABdA=zeros(nAB,nA);
end;

ni = 1:p:p*(q-1)+1;
nj = 1:p:p*(n-1)+1;
for i=1:p,
    dABdA(ni+i-1,nj+i-1) = B';
end;

if issparse(A) ||  issparse(B) || nAB*nB>400,
  dABdB = speye(q);
else
  dABdB = eye(q);
end;
dABdB = kron(dABdB,A);        % Kronecker tensor product.

return;


%% other methods
% for j=1:q,
%    for i=1:p,
%    ij = i + (j-1)*p;
%       for k=1:n,
%          ik = i + (k-1)*p;
%          dABdA(ij,ik) = B(k,j);
%       end;
%    end;
% end;

% if nargout <2,
%     return;
% end;
% 
% if issparse(A) ||  issparse(B) || nAB*nB>400,
%   dABdB=sparse([],[],[],nAB,nB,q*nA);
% else
%   dABdB=zeros(nAB,nB);
% end;
% for j=1:q
%    dABdB([j*p-p+1:j*p]',[j*n-n+1:j*n]) = A;
% end;




%%  Test
p1 = 8;
q1 = 7;
n1 = 6;
a = randn(p1,n1);
b = randn(n1,q1);
ab = a*b;
[dabda, dabdb] = dAB(a, b);
da = randn(p1,n1)/1000;
db = randn(n1,q1)/1000;
ab1 = (a+da)*(b+db);
ab2 = ab+reshape(dabda*da(:),p1,q1)+reshape(dabdb*db(:),p1,q1);
gain = norm(ab1-ab)/norm(ab1-ab2)

