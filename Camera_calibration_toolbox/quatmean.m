function [ Qm ] = quatmean( Qs, Ws )
% QUATMEAN - The mean quaternion of a bunch of weighted unit quaternions. Qs are
%   quaternions with 4*n elements. 1~3 rows are imaginary parts (3d vector part),
%   the 4th row is real part (scalor part). Ws are the weights corresponding to
%   input quaternions.
%
%   reference: Markley, Averaging Quaternions, Journal of guidance, control, and
%   dynamics, 2007, p1194.

%  by zpf, form BIT, 2016-12-22

[m,n] = size(Qs);
assert(ismatrix(Qs) && m==4,'The 1st argument must be a matrix with 4 rows!');
assert(isvector(Ws) && length(Ws)==n && all(Ws>=0),['The 2nd argument must be '...
      'a non-negative vector with same columns with 1st argument!']);

Qs2 = Qs.*Qs;
T = sqrt(sum(Qs2));
Qs = Qs./T(ones(4,1),:);
Qs2 = zeros(4,4,n);
for i=1:n,
  Qs2(:,:,i) = Ws(i)*(Qs(:,i)*Qs(:,i)');
end;
Qs2 = sum(Qs2,3);
[~,~,T] = svd(Qs2);
Qm = T(:,1);

return;


%% Test
n = 2;
a = randn(4,n);
a = a./(ones(4,1)*sqrt(sum(a.^2)));
b = ones(1,n)/n;
%b = rand(1,n); b = b/sum(b);
q = quatmean(a,b);
for i=2:n,
  if sum(prod(a(:,[1,i]),2))<0,
    a(:,i)=-a(:,i);
  end;
end;
c=sum(a,2)/n;
c=c/norm(c);
