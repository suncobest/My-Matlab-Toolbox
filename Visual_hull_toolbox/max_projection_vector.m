function [d, vec, ind] = max_projection_vector(X,V)
% MAX_PROJECTION_VECTOR - compute maximum projection (on vector V) vector
% composed by two points in X.
%
%    X - (DxM) points matrix; 
%    V - (Dx1) given vector;
% 
% d: projection length on V of vector vec;
% vec: maximum projection vector;
% ind: index of two ends in data X;
% 
% [d,vec,ind] = max_projection_vector(X) is the same as
% [d,vec,ind] = max_distance(X,X);
%
% Example : 
%    X = rand(3,20); V = randn(3,1);
%    [d, vec, ind] = max_projection_vector(X,V);
% See also max_distance, min_distance, L2_distance.

% by zpf, from BIT, Jan 20, 2016.

assert(isreal(X),'Points must be real!'); 
[nd, np] = size(X);
if nargin < 2 || isempty(V),
    aa = repmat(sum(X.^2,1),[np,1]);
    E = sqrt(aa'+aa- 2*(X'*X));
    [d, I] = max(E,[],1);
    [d, J] = max(d,[],2);
    I = I(J);
    ind = [I,J];
    vec = X(:,J)-X(:,I);
else
    assert(isreal(V) && isvector(V) && length(V)==nd, 'Unexpected input for the 2nd argument!');
    X = X-repmat(mean(X,2),1,np);
    V = V(:)/norm(V(:));
    Y = V'*X;
    [a,I] = min(Y);
    [b,J] = max(Y);
    d = b-a;
    ind = [I,J];
    vec = X(:,J)-X(:,I);
end;

return;


%% Test

a = 10*randn(2,30);
b = randn(2,1);
c = mean(a,2);
cb = [c, c+b/norm(b)*20];
[d, vc, id] = max_projection_vector(a,b);
err = d-vc'*b/norm(b)
figure;
a12 = [a(:,id(1)), a(:,id(2))];
plot(a(1,:),a(2,:),'.', c(1),c(2),'+', cb(1,:),cb(2,:),'k-', a(1,id(1)),a(2,id(1)),'o', a(1,id(2)),a(2,id(2)),'o', a12(1,:), a12(2,:),'r-');
axis image off; hold on;
[d, vc, id] = max_projection_vector(a);
a12 = [a(:,id(1)), a(:,id(2))];
plot(a(1,id(1)),a(2,id(1)),'o', a(1,id(2)),a(2,id(2)),'o', a12(1,:), a12(2,:),'g-');

a = 10*randn(3,100);
b = randn(3,1);
c = mean(a,2);
cb = [c, c+b/norm(b)*20];
[d, vc, id] = max_projection_vector(a,b);
err = d-vc'*b/norm(b)
figure;
a12 = [a(:,id(1)), a(:,id(2))];
plot3(a(1,:),a(2,:),a(3,:),'.', c(1),c(2),c(3),'+', cb(1,:),cb(2,:),cb(3,:),'k-', a(1,id(1)),a(2,id(1)),a(3,id(1)),'o', a(1,id(2)),a(2,id(2)),a(3,id(2)),'o', a12(1,:), a12(2,:),a12(3,:),'r-');
axis image off; hold on;
[d, vc, id] = max_projection_vector(a);
a12 = [a(:,id(1)), a(:,id(2))];
plot3(a(1,id(1)),a(2,id(1)),a(3,id(1)),'o', a(1,id(2)),a(2,id(2)),a(3,id(2)),'o', a12(1,:), a12(2,:),a12(3,:),'g-');

