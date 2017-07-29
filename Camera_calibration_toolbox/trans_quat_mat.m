function [out, dout] = trans_quat_mat( in )
% TRANS_QUAT_MAT - transform quaternion to matrix, or vice versa.
%  QUAT2MAT transform quaternion to matrix, while MAT2QUAT transform matrix to quaternion.
%
%  Assuming that a quaternion has the form:
%
%    Q = |X Y Z W|
%
%   norm([X Y Z]) = sin(theta/2), W = cos(theta/2)
%
%  Then the quaternion can be converted into a 3x3 rotation matrix using the following expression:
%
%         | 1-2Y^2-2Z^2     2XY - 2ZW        2XZ + 2YW   |
%    R =  | 2XY + 2ZW       1-2X^2-2Z^2      2YZ - 2XW   |
%         | 2XZ - 2YW       2YZ + 2XW        1-2X^2-2Y^2 |
%
%   dout is the derivative of the output wrt the input
%   wrt : with respect to
%
%   Algorithm:
%   Advanced Animation And Rendering Techniques - Theory And Practice, P362-365
%   See also trans_quat_axis, rodrgues, slerp.

% by zpf, from BIT, 2015-7-12

bigeps = 1e5;     %1e+16*eps;

if isvector(in) && length(in)==4,
    if nargout < 2,
        out = quat2mat(in);
    else
        [out, dout] = quat2mat(in);
    end;
elseif  isequal(size(in),[3,3]) && norm(in' * in - eye(3)) < bigeps,
    detR = det(in);
    assert( detR>0 && abs(detR-1)<bigeps, 'The input rotation matrix (3x3) must have positive determinant!');
    if nargout < 2,
        out = mat2quat(in);
    else
        [out, dout] = mat2quat(in);
    end;
else
    error('The input argument must either be a quaternion or a rotation matrix!');
end;

return;


function [R, dRdQ] = quat2mat(Q)

T = norm(Q);
assert(T>0, 'The quaternion can not be a zero vector!');
Q=Q/T;    % ensure unit quaternion
X=Q(1);
Y=Q(2);
Z=Q(3);
W=Q(4);

X2 = X^2;
Y2 = Y^2;
Z2 = Z^2;
W2 = W^2;
XY = X*Y;
XZ = X*Z;
XW = X*W;
YZ = Y*Z;
YW = Y*W;
ZW = Z*W;

R = [1-2*(Y2+Z2),  2*(XY-ZW),  2*(XZ+YW);
    2*(XY+ZW),  1-2*(X2+Z2),  2*(YZ-XW);
    2*(XZ-YW),  2*(YZ+XW),  1-2*(X2+Y2)];

if nargout < 2,
    return;
end;

% dR/dQ1
dRdQ = [0, -4*Y, -4*Z, 0;
    2*Y, 2*X, 2*W, 2*Z;
    2*Z, -2*W, 2*X, -2*Y;
    2*Y, 2*X, -2*W, -2*Z;
    -4*X, 0, -4*Z, 0;
    2*W, 2*Z, 2*Y, 2*X;
    2*Z, 2*W, 2*X, 2*Y;
    -2*W, 2*Z, 2*Y, -2*X;
    -4*X, -4*Y, 0, 0];

% dQ1/dQ0, det(dQ)=0
% dQ = [1-X2, -XY, -XZ, -XW;
%     -XY, 1-Y2, -YZ, -YW;
%     -XZ, -YZ, 1-Z2, -ZW;
%     -XW, -YW, -ZW, 1-W2]/T;
% dRdQ = dR*dQ;
return;


function [Q, dQdR] =mat2quat(R)

[U,~,V] = svd(R);
R = U*V';       % The closest rotation matrix

T = sum(diag(R));       % 4W^2-1
Q = zeros(4,1);

if ( T > 0 ),
    S = 0.5 / sqrt(T+1);
    Q(4) = 0.25 / S;
    Q1 = R(3,2) - R(2,3);
    Q2 = R(1,3) - R(3,1);
    Q3 = R(2,1) - R(1,2);
    Q(1:3) = [Q1; Q2; Q3] * S;
    if nargout >1,
        dTdR = [1, 0, 0, 0, 1, 0, 0, 0, 1];
        dSdT = -2*S^3;
        dSdR = dSdT*dTdR;
        dQ1dR = [0, 0, 0, 0, 0, 1, 0, -1, 0];
        dQ2dR = [0, 0, -1, 0, 0, 0, 1, 0, 0];
        dQ3dR = [0, 1, 0, -1, 0, 0, 0, 0, 0];
        dQdR = zeros(4, 9);
        dQdR(4,:) = 0.5*S*dTdR;
        dQdR(1:3, :) = [dQ1dR; dQ2dR; dQ3dR]*S + [Q1; Q2; Q3] * dSdR;
    end;
else                          % abs(W)<=1/2
    [~, i] = max(diag(R));
    nxt = [2 3 1];
    j = nxt(i);
    k = nxt(j);
    S  = 0.5/sqrt( 1.0 + R(i,i) - R(j,j) - R(k,k) );
    Q(i) = 0.25 / S;
    Qj = R(i,j) + R(j,i);
    Qk = R(i,k) + R(k,i);
    Q4 = R(k,j) - R(j,k);
    Q(j) = Qj*S;
    Q(k) = Qk*S;
    Q(4) = Q4*S;
    if nargout >1,
        dTdR = zeros(1,9);
        dQjdR = dTdR;
        dQkdR = dTdR;
        dQ4dR = dTdR;
        dTdR(4*i-3) = 1;
        dTdR(4*j-3) = -1;
        dTdR(4*k-3) = -1;
        dSdT = -2*S^3;
        dSdR = dSdT*dTdR;
        dQdR = zeros(4, 9);
        dQdR(i,:) = 0.5*S*dTdR;
        dQjdR(3*(j-1)+i) =1;
        dQjdR(3*(i-1)+j) =1;
        dQkdR(3*(k-1)+i) =1;
        dQkdR(3*(i-1)+k) =1;
        dQ4dR(3*(j-1)+k) =1;
        dQ4dR(3*(k-1)+j) =-1;
        dQdR(j,:) = dQjdR*S + Qj * dSdR;
        dQdR(k,:) = dQkdR*S + Qk * dSdR;
        dQdR(4,:) = dQ4dR*S + Q4 * dSdR;
    end;
end;

return;


%% Test
q = randn(4,1);
q = q/norm(q);
q1=q+randn(4,1)/300;
dq = q1/norm(q1)-q;
[r, drdq] =trans_quat_mat(q);  % drdq is the derivative of r wrt q(矩阵对向量求导)
r1 = trans_quat_mat(q+dq);
r2 = r + reshape(drdq * dq,3,3);
gain = norm(r1 - r)/norm(r1- r2) ;  % 分子是函数变化量，分母是函数余项（二阶以上的高阶项）。大量/小量

%% Test of drdq:
gain = 1000;
i=0;
while gain>100,
    i=i+1;
    q = randn(4,1);
    q = q/norm(q);
    q1=q+randn(4,1)/300;
    dq = q1/norm(q1)-q;
    [r, drdq] =trans_quat_mat(q);  % drdq is the derivative of r wrt q(矩阵对向量求导)
    r1 = trans_quat_mat(q+dq);
    r2 = r + reshape(drdq * dq,3,3);
    gain = norm(r1 - r)/norm(r1- r2) ;  % 分子是函数变化量，分母是函数余项（二阶以上的高阶项）。大量/小量
end;

%% Test of dqdr:
gain = 1000;
i=0;
while gain>100,
    i=i+1;
    om = randn(3,1);
    r = rodrigues(om);
    dom = randn(3,1)/100;
    dr = rodrigues(om+dom)-r;
    [q, dqdr] =trans_quat_mat(r);
    q1 = trans_quat_mat(r+dr);
    if sum(q1.*q)<0,
        q1=-q1;
    end;
    q2 = q + dqdr * dr(:);
    q2 = q2/norm(q2);
    gain = norm(q1 - q)/norm(q1- q2);   % 大量/小量
end;
