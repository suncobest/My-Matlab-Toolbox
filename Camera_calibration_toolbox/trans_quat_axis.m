function [ out, dout ] = trans_quat_axis( in )
% TRANS_QUAT_AXIS - transform quaternion to axis-angle, or vice versa.
%  QUAT2AXIS transform quaternion to axis-angle, while AXIS2QUAT transform axis-angle to quaternion.
%
%  Assuming that a quaternion has the form:
%
%    Q = |X Y Z W|
%
%   norm([X Y Z]) = sin(theta/2), W = cos(theta/2)
%
%  Then the corresponding axis-angle can be  expressed in the form:
%
%  OM = theta*[X Y Z]/norm([X Y Z])
%
%  dout is the derivative of the output wrt the input
%  wrt : with respect to
%   See also trans_quat_mat, rodrgues, slerp.

% by zpf, form BIT, 2015-7-12

assert(ismatrix(in),'Input must be a vector or matrix!');
if isvector(in),
    switch length(in),
        case 4,
            in = in(:);
            if nargout < 2,
                out = quat2axis(in);
            else
                [out,dout] = quat2axis(in);
            end;
        case 3,
            in = in(:);
            if nargout < 2,
                out = axis2quat(in);
            else
                [out, dout] = axis2quat(in);
            end;
        otherwise,
            error('Input arguments must be either quaternions or axis-angles!');
    end;
else
    ndim = size(in,1);
    switch ndim,
        case 4,
            if nargout < 2,
                out = quat2axis(in);
            else
                [out,dout] = quat2axis(in);
            end;
        case 3,
            if nargout < 2,
                out = axis2quat(in);
            else
                [out, dout] = axis2quat(in);
            end;
        otherwise,
            error('Input arguments must be either quaternions or axis-angles!');
    end;
end;

return;


function [OM, dOMdQ] = quat2axis(Q)

% make sure W>=0, so 0<=theta<=pi
ids = sign(Q(4,:))<0;
Q(:,ids) = -Q(:,ids);
T = sqrt(sum(Q.*Q,1));
assert(all(T>0), 'Quaternions can not be zero vectors!');
Q1 = Q./T(ones(4,1),:);    % ensure unit quaternion
W = Q1(4,:);
one_W2 = 1-W.^2;

theta = acos(W)*2;              % 0<=theta<=pi
stheta = sqrt(one_W2);       % sin(theta/2)
OM1 = Q1(1:3,:);
n1 = ones(3,1);
id = stheta>=1e-7;           % sin(theta/2)>=1e-7, so 2e-7<=theta<=pi
OM = OM1*2;       % if 0<=theta<2e-7, sin(theta/2) = theta/2
OM(:,id) = theta(n1,id).*OM1(:,id)./stheta(n1,id);

if nargout < 2,
    return;
end;

X=Q1(1,:);
Y=Q1(2,:);
Z=Q1(3,:);
one_X2 = 1-X.^2;
one_Y2 = 1-Y.^2;
one_Z2 = 1-Z.^2;
n = length(X);
dOMdQ = zeros(3,4,n);
XY = X.*Y;
XZ = X.*Z;
XW = X.*W;
YZ = Y.*Z;
YW = Y.*W;
ZW = Z.*W;
for i = 1:n,
    dQ1dQ = [one_X2(i), -XY(i), -XZ(i), -XW(i);
        -XY(i), one_Y2(i), -YZ(i), -YW(i);
        -XZ(i), -YZ(i), one_Z2(i), -ZW(i);
        -XW(i), -YW(i), -ZW(i), one_W2(i)]/T(i);
    if ids(i),
        dQ1dQ = -dQ1dQ;
    end;
    % dOMdQ1    :  Q1=[OM1; W]
    if id(i),     % OM = OM1*2*acos(W)/sqrt(1-W^2)
        dOMdOM1 = eye(3)*theta(i)/stheta(i);
        dOMdW = (OM(:, i)*W(i)-2*OM1(:, i))/one_W2(i);
    else
        dOMdOM1 = eye(3)*2;
        dOMdW =  [0; 0; 0];
    end;
    dOMdQ(:,:,i) = [dOMdOM1, dOMdW]*dQ1dQ;
end;

return;


function [Q, dQdOM] = axis2quat(OM)

theta = sqrt(sum(OM.*OM, 1));
W = cos(theta/2);
stheta = sin(theta/2);
n1 = ones(3,1);

id = theta>1e-7;
OM1 = OM/2;      % if 0<=theta<=1e-7,  sin(theta/2) = theta/2
OM1(:,id) = stheta(n1,id).*OM(:,id)./theta(n1,id);
Q = [OM1; W];

if nargout < 2,
    return;
end;
n = length(W);
dQdOM = zeros(4,3,n);
for i=1:n,
    if id(i),
        dthetadOM = OM(:,i)'/theta(i);
        dWdtheta = -0.5*stheta(i);
        dWdOM = dWdtheta*dthetadOM;
        dOM1dOM = eye(3)*stheta(i)/theta(i)+(0.5*W(i)*OM(:,i)-OM1(:,i))/theta(i)*dthetadOM;
    else
        dWdOM = -OM(:,i)'/4;
        dOM1dOM = 0.5*eye(3);
    end;
    dQdOM(:,:,i) = [dOM1dOM; dWdOM];
end;

return;


%% Test of domdq:
gain = 1000;
i=0;
while gain >100,
    i=i+1;
    q = randn(4,1);
    dq = randn(4,1)/1000;
    [om, domdq] =trans_quat_axis(q);  % domdq is the derivative of om wrt q
    om1 = trans_quat_axis(q+dq);
    if norm(om1-om)>1,
        t1 = norm(om1);
        om1 = om1/t1;
        om1 = (t1-2*pi)*om1;
    end;
    om2 = om + domdq * dq;
    gain = norm(om1 - om)/norm(om1- om2);   % O(dq)/O(dq^2)
end;

%% Test of dqdom:
gain = 1000;
i=0;
while gain >100,
    i=i+1;
    om = randn(3,1);
    dom = randn(3,1)/100;
    [q, dqdom] =trans_quat_axis(om);
    q1 = trans_quat_axis(om+dom);
    q2 = q + dqdom * dom;
    gain = norm(q1 - q)/norm(q1- q2);   % big/small
end;
