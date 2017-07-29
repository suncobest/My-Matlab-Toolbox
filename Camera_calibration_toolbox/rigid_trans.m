function [Y,dYdQ,dYdT,dYdX] = rigid_trans(X,Q,T,hand)
% rotation and refmotion
% ref: rigid_refmotion.m
%
% [Y,dYdX,dYdom,dYdT] = rigid_trans(X,om,T,hand)
% Computes the rigid transformation: Y = quatmul(quatmul(Q, diag([1 1 hand]*X), invQ)+T;
% The algorithm is equivalent to: Y= R*X+T, where R = trans_quat_mat(Q)*diag([1 1 hand]).
%
% INPUT: X: 3D structure in the world coordinate frame (3xN or 2xN matrix for N points)
%       (Q,T): Rigid transformation parameters between world coordinate frame and camera reference frame
%       Q: rotation quaternion (4x1 vector); T: translation vector (3x1vector); hand: handness of frame (1 or -1)
%
% OUTPUT: Y: 3D coordinates of the structure points in the camera reference frame (3xN matrix for N points)
%        dYdQ: Derivative of Y with respect to Q ((3N)x4 matrix)
%        dYdT: Derivative of Y with respect to T ((3N)x3 matrix)
%
% Definitions:
% Let P be a point in 3D of coordinates X in the world reference frame (stored in the matrix X)
% The coordinate vector of P in the camera reference frame is: Y = R*X + T
% where R is the orthogonal matrix corresponding to the rotation vector om:
% R = trans_quat_mat(Q)*diag([1 1 hand]);
% See also rigid_refmotion, rigid_motion, quatmul, trans_quat_mat, trans_quat_axis, rodrgues.

if nargin < 4,
    hand = 1;
    if nargin < 3,
        T = zeros(3,1);
        if nargin < 2,
            Y = X;
            dYdQ = [];
            dYdT = [];
            dYdX = speye(numel(X));
            return;
        end;
    end;
end;
[m,n] = size(T);
assert(m==3 && n==1, 'The 3rd argument (translation) must have 3 rows and 1 column!');
[m,n] = size(X);
if m==2,
    X = [X;zeros(1,n)];
else
    assert(m==3, 'Unexpected dimension for the 1st argument!');
end;
nQ = norm(Q);
assert(nQ>0, 'The 2nd argument (quaternion) can not be a zero vector!');
Q=Q/nQ;

q1=Q(1);
q2=Q(2);
q3=Q(3);
q4=Q(4);

X2 = q1^2;
Y2 = q2^2;
Z2 = q3^2;
% W2 = q4^2;
XY = q1*q2;
XZ = q1*q3;
XW = q1*q4;
YZ = q2*q3;
YW = q2*q4;
ZW = q3*q4;

R = [1-2*(Y2+Z2),  2*(XY-ZW),  2*(XZ+YW);
    2*(XY+ZW),  1-2*(X2+Z2),  2*(YZ-XW);
    2*(XZ-YW),  2*(YZ+XW),  1-2*(X2+Y2)];

if hand ~= 1,
    R(:,3)  = -R(:,3);
end;
Y = R*X + repmat(T,[1 n]);

if nargout < 2,
    return;
end;

% dR/dQ1
dRdQ = [0, -4*q2, -4*q3, 0;
    2*q2, 2*q1, 2*q4, 2*q3;
    2*q3, -2*q4, 2*q1, -2*q2;
    2*q2, 2*q1, -2*q4, -2*q3;
    -4*q1, 0, -4*q3, 0;
    2*q4, 2*q3, 2*q2, 2*q1;
    2*q3, 2*q4, 2*q1, 2*q2;
    -2*q4, 2*q3, 2*q2, -2*q1;
    -4*q1, -4*q2, 0, 0];

% dQ1/dQ0, det(dQ)=0
% dQ = [1-X2, -XY, -XZ, -XW;
%     -XY, 1-Y2, -YZ, -YW;
%     -XZ, -YZ, 1-Z2, -ZW;
%     -XW, -YW, -ZW, 1-W2]/nQ;
% dRdQ = dRdQ*dQ;

if hand ~= 1,
	dRdQ(7:9,:) = -dRdQ(7:9,:);
end;
n3 = n*3;
if n3*3>450;
    dYdR = sparse([],[],[],n3,9,n3*3);
else
    dYdR = zeros(n3,9);
end;

dYdR(1:3:end,1:3:end) =  X';
dYdR(2:3:end,2:3:end) =  X';
dYdR(3:3:end,3:3:end) =  X';
if m==2,
    dRdQ = dRdQ(1:6,:);
    dYdR = dYdR(:,1:6);
end;
dYdQ = dYdR * dRdQ;

dYdT = zeros(n3,3);
dYdT(1:3:end,1) =  ones(n,1);
dYdT(2:3:end,2) =  ones(n,1);
dYdT(3:3:end,3) =  ones(n,1);

if nargout < 4,
    return;
end;

if n>5;
    dYdX = speye(n);
else
    dYdX = eye(n);
end;
if m==2,
    dYdX = kron(dYdX,R(:,1:2));
else
    dYdX = kron(dYdX,R);        % Kronecker tensor product.
end;

return;



%% test
m = 2;
om=randn(3,1);
om=om/norm(om);
a=pi-rand/100000;
q=[om*sin(a/2);cos(a/2)];
om = om*a;
np = 100;
nd = 100;
h = -1;
x=randn(m,np)*10;
t = randn(3,1)*10;

dx = randn(m,np)/nd;
dom = randn(3,1)/nd;
dt = randn(3,1)/nd;    % zeros(3,1);
dq = trans_quat_axis(om+dom)-q;
% q1 = q+randn(4,1)/nd;
% dq=q1/norm(q1)-q;

[y, dydom,dydt,dydx] = rigid_refmotion(x,om,t,h);
y1 = rigid_refmotion(x+dx,om+dom,t+dt,h);
y2 = y +reshape(dydx*dx(:),3,[])+ reshape(dydom*dom(:),3,[]) + reshape(dydt*dt(:),3,[]);
gain = norm(y1-y)/norm(y1-y2)

[y, dydq,dydt,dydx] = rigid_trans(x,q,t,h);
y1 = rigid_trans(x+dx,q+dq,t+dt,h);
y2 = y + reshape(dydx*dx(:),3,[]) + reshape(dydq*dq(:),3,[]) + reshape(dydt*dt(:),3,[]);
gain = norm(y1-y)/norm(y1-y2)

