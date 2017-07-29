function [d, vec, ind, E] = min_distance(a,b)
% min_DISTANCE - computes mininum L2 Euclidean distance of two points cloud.
%
% D = min_distance(A,B)
%
%    A - (DxM) matrix 
%    B - (DxN) matrix
% 
% Returns: Min Euclidean distance between vectors in A and B
%
% Description : 
%    This fully vectorized (VERY FAST!) m-file computes the 
%    Euclidean distance between two vectors by:
%
%                 ||A-B|| = sqrt ( ||A||^2 + ||B||^2 - 2*A.B )
%
% Example : 
%    A = rand(400,100); B = rand(400,200);
%    d = min_distance(A,B);
% See also max_distance, L2_distance, max_projection_vector.

% Author   : Roland Bunschoten
%            University of Amsterdam
%            Intelligent Autonomous Systems (IAS) group
%            Kruislaan 403  1098 SJ Amsterdam
%            tel.(+31)20-5257524
%            bunschot@wins.uva.nl
% Last Rev : Wed Oct 20 08:58:08 MET DST 1999
% Tested   : PC Matlab v5.2 and Solaris Matlab v5.3

% Copyright notice: You are free to modify, extend and distribute 
%    this code granted that the author of the original code is 
%    mentioned as the original author of the code.
% Modified by zpf, from BIT, Jan 20, 2016.

if nargin < 2,
    error('Not enough input arguments');
else
    assert(size(a,1) == size(b,1), 'A and B should be of same dimensionality');
end;
assert(isreal(a) && isreal(b),'Points must be real!'); 

aa=sum(a.*a,1); bb=sum(b.*b,1); ab=a'*b; 
E = sqrt(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);
[d, I] = min(E,[],1);
[d, J] = min(d,[],2);
I = I(J);
ind = [I,J];
vec = b(:,J)-a(:,I);
return;


%% Test
a = randn(2,200);
b = randn(2,300)+5;
[d, vc, id] = min_distance(a,b);
figure;
lab = [a(:,id(1)), b(:,id(2))]';
plot(a(1,:),a(2,:),'.', b(1,:),b(2,:),'.', a(1,id(1)),a(2,id(1)),'o', b(1,id(2)),b(2,id(2)),'o', lab(:,1), lab(:,2),'-');
axis image off;
err = d-norm(vc)

a = 10*randn(3,500);
b = 10*randn(3,300)+50;
[d, vc, id] = min_distance(a,b);
figure;
lab = [a(:,id(1)), b(:,id(2))]';
plot3(a(1,:),a(2,:),a(3,:),'.', b(1,:),b(2,:),b(3,:),'.', a(1,id(1)),a(2,id(1)),a(3,id(1)),'o', b(1,id(2)),b(2,id(2)),b(3,id(2)),'o', lab(:,1), lab(:,2),lab(:,3),'-');
axis image off;
err = d-norm(vc)

