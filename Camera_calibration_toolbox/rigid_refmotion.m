function [Y,dYdom,dYdT,dYdX] = rigid_refmotion(X,om,T,hand)
%ref:    rigid_motion.m
%
%[Y,dYdom,dYdT] = rigid_refmotion(X,om,T,hand)
%
%Computes the rigid motion transformation Y = R*X+T, where R = rodrigues(om)*diag([1 1 hand]).
%
%INPUT: X: 3D structure in the world coordinate frame (3xN or 2xN matrix for N points)
%       (om,T): Rigid motion parameters between world coordinate frame and camera reference frame
%               om: rotation vector (3x1 vector); T: translation vector (3x1 vector)
%
%OUTPUT: Y: 3D coordinates of the structure points in the camera reference frame (3xN matrix for N points)
%        dYdom: Derivative of Y with respect to om ((3N)x3 matrix)
%        dYdT: Derivative of Y with respect to T ((3N)x3 matrix)
%
%Definitions:
%Let P be a point in 3D of coordinates X in the world reference frame (stored in the matrix X)
%The coordinate vector of P in the camera reference frame is: Y = R*X + T
%where R is the rotation matrix corresponding to the rotation vector om: R = rodrigues(om);
%
%Important function called within that program:
%
%rodrigues.m: Computes the rotation matrix corresponding to a rotation vector

if nargin < 4,
    hand = 1;
    if nargin < 3,
        T = zeros(3,1);
        if nargin < 2,
            Y = X;
            dYdom = [];
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

[R,dRdom] = rodrigues(om);
if hand ~= 1,
    R(:,3)  = -R(:,3);
    dRdom(7:9,:) = -dRdom(7:9,:);
end;

Y = R*X + repmat(T,[1 n]);

if nargout < 2,
    return;
end;
n3 = n*3;
dYdR = zeros(n3,9);

dYdR(1:3:end,1:3:end) =  X';
dYdR(2:3:end,2:3:end) =  X';
dYdR(3:3:end,3:3:end) =  X';
if m==2,
    dRdom = dRdom(1:6,:);
    dYdR = dYdR(:,1:6);
    R=R(:,1:2);
end;
dYdom = dYdR * dRdom;

dYdT = zeros(n3,3);
dYdT(1:3:end,1) =  ones(n,1);
dYdT(2:3:end,2) =  ones(n,1);
dYdT(3:3:end,3) =  ones(n,1);

if nargout < 4,
    return;
end;

if n>5;
    dYdX = sparse([],[],[],n3,n*m,n3*m);
else
    dYdX = zeros(n3,n*m);
end;

for i=1:n,
    ii = (i-1)*3;
    jj = (i-1)*m;
    dYdX(ii+1:ii+3,jj+1:jj+m) = R;
end;

return;



%% Test of jacobian
m1 = randi(2)+1
n1 = 1000;
X1 = randn(m1,n1);
om1 = randn(3,1);
T1 = randn(3,1);
h1 = sign(randn);
[Y1,dYdom1,dYdT1,dYdX1] = rigid_refmotion(X1,om1,T1,h1);
dX1 = randn(m1,n1)/100;
dom1 = randn(3,1)/100;
dT1 = randn(3,1)/100;
Y2 = rigid_refmotion(X1+dX1,om1+dom1,T1+dT1,h1);
Yp = Y1+reshape(dYdX1*dX1(:)+dYdom1*dom1+dYdT1*dT1,3,n1);
gain = norm(Y2-Y1)/norm(Y2-Yp)

% dYdX
Y2 = rigid_refmotion(X1+dX1,om1,T1,h1);
Yp = Y1+reshape(dYdX1*dX1(:),3,n1);
gain = norm(Y2-Y1)/norm(Y2-Yp)

% dYdom
Y2 = rigid_refmotion(X1,om1+dom1,T1,h1);
Yp = Y1+reshape(dYdom1*dom1,3,n1);
gain = norm(Y2-Y1)/norm(Y2-Yp)

% dYdT
Y2 = rigid_refmotion(X1,om1,T1+dT1,h1);
Yp = Y1+reshape(dYdT1*dT1,3,n1);
gain = norm(Y2-Y1)/norm(Y2-Yp)
