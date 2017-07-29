function [index, Xout] = points_between_2planes( Xin, P1, N1, P2, N2 )
%POINTS_BETWEEN_2PLANES compute points between two planes.
%   [index, Xout] = points_between_2planes( Xin, P1, N1, P2, N2 )
%   Compute points between 2 planes. (>=plane1 and <=plane2)
%   Xin: input 3D points of dimension (3, n);
%   Xout: output 3D points of dimension (3, n);
%   P1: a point on plane 1;
%   P2: a point on plane 2;
%   N1: normal of plane 1;
%   N2: normal of plane 2;
%   index: logical column index of 3D points;
%
%   [index, Xout] = points_between_2planes( Xin, P1, N1 )
%   compute points on one side of plane1 (>=);

[m, n] = size(Xin);
assert(m==3, 'Unexpected dimension of input points!');

if nargin<3,
    N1 = [0; 0; 1];
    if nargin<2,
        P1 = [0; 0; 0];
    end;
end;
assert(isequal(size(P1),[3,1]),'Unexpected dimension of point on plane1!');
assert(isequal(size(N1),[3,1]),'Unexpected dimension of plane 1 normal!');
index = sum((Xin-P1(:,ones(1,n))).*N1(:,ones(1,n)))>=0;
Xout = Xin(:,index);

if nargin <4,
    return;
end;
if nargin<5,
    N2 = [0; 0; 1];
end;
assert(isequal(size(P2),[3,1]),'Unexpected dimension of point on plane2!');
assert(isequal(size(N2),[3,1]),'Unexpected dimension of plane 2 normal!');
index = sum((Xin-P2(:,ones(1,n))).*N2(:,ones(1,n)))<=0 & index;
Xout = Xin(:,index);

return;


%% Test
XX = 10*diag([3 2 1])*randn(3,10000);
p1 = 5*randn(3,1);
p2 = 5*randn(3,1);
n1 = randn(3,1);
n2 = randn(3,1);
n1 = n1/norm(n1);
n2 = n2/norm(n2);
id1 = points_between_2planes(XX,p1,n1,p2,n2);
id2 = points_between_2planes(XX,p1,n1);
figure;
plot3(XX(1,id1),XX(2,id1),XX(3,id1),'.',XX(1,~id1),XX(2,~id1),XX(3,~id1),'.');
axis image off;
figure;
plot3(XX(1,id2),XX(2,id2),XX(3,id2),'.',XX(1,~id2),XX(2,~id2),XX(3,~id2),'.');
axis image off;
