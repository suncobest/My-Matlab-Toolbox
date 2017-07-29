function [om,T,hand] = compute_rigid_refmotion(Y,X,SW)

% COMPUTE_RIGID_REFMOTION computes the rigid motion (rotation, translation and reflection).
%
%        Y = R*X+T, where R = rodrigues(om)*diag([1 1 hand])
%
% INPUT: X: The 3D coordinates of rigid structure points under body-fixed frame of axes; (3xN)
%        Y: The corresponding points measured with Gaussian white noise after motion; (3xN)
%        SW: switch turning on/off the optimization process.
%
% OUTPUT: (om,T,hand): Motion parameters of the 3D rigid structure;
%        om: rotation vector (3x1 vector);
%        T: translation vector (3x1 vector);
%        hand: the handedness of the orthogonal matrix R; R is rotation if hand==1;
%              R is rotation after reflection if hand==-1;
%
% ALGORITHM: The position of geometric center is more accurate than individual points in general.
%        So we use the geometric center as reference to compute the initial value of translation.
%        The furthest point from centroid is used to build axes to locate the body frame's origin.
%
%        The rotation is then computed using SVD decomposition.
%
%  See also compute_rigid_refmotion2, rigid_refmotion, compose_motion2.

% By ZPF @ZVR, 2017-7-21

if nargin<3, SW=1; end;
[my,ny] = size(Y);
[m,n] = size(X);
assert(m==3 && my==3, 'The input matrix must have 3 rows!');
assert(n==ny, 'The two input matrix must have the same number of points!');
assert(n>=3, 'It takes at least 3 points to determine a rigid body!');
bigeps = eps*1e6;

% Schmidt orthogonalization
Cx = mean(X,2);
XXc = X-Cx(:,ones(1,n));
[nax2,ia] = max(sum(XXc.^2,1));
Cax = XXc(:,ia)/sqrt(nax2);
rem = XXc-Cax*(Cax'*XXc);
ind = 1:n;
ind(ia) = [];
[nbx2,i] = max(sum(rem(:,ind).^2,1));
ib = ind(i);
ind(i) = [];
Cbx = rem(:,ib)/sqrt(nbx2);
Rx = [Cax,Cbx,cross(Cax,Cbx)];
remz = Rx(:,3)'*rem;
[~,i] = max(abs(remz(ind)));
ic = ind(i);

if remz(ic)>bigeps,
    zsx = 1;
elseif remz(ic)<-bigeps,
    zsx = -1;
else
    zsx = 0;
end;

% coordinates of body frame's origin in Cartesian system (Cax,Cbx,cross(Cax,Cbx))
lamda = -Rx'*Cx;

Cy = mean(Y,2);
YYc = Y(:,[ia,ib,ic])-Cy(:,ones(1,3));
Cay = YYc(:,1)/norm(YYc(:,1));
Cby = YYc(:,2)-Cay*(Cay'*YYc(:,2));
Cby = Cby/norm(Cby);
Ry = [Cay,Cby,cross(Cay,Cby)];

remz = Ry(:,3)'*YYc(:,3);
if remz>bigeps,
    zsy = 1;
elseif remz<-bigeps,
    zsy = -1;
else
    zsy = 0;
end;

if zsy == zsx,
    hand = 1;
else
    hand = -1;
    Ry(:,3) = -Ry(:,3);
end;
% another method: compute rotation before translation
% R = Ry*Rx';
% om = rodrigues(R*diag([1,1,hand]));
% T = mean(Y-R*X,2);

% Initial value of translation
T = Cy+Ry*lamda;
[u,~,v] = svd((Y-T(:,ones(1,n)))*X');

% Initial value of rotation
R = u*v';
if hand == -1,
    R(:,3) = -R(:,3);
end;
om = rodrigues(R);

if ~SW, return; end;

% Optimization through LM algorithm
param = [om;T];
lamda = 0.001; % set an initial value of the damping factor for the LM
updateJ = 1;
Y_up = rigid_refmotion(X,om,T,hand);
err = Y - Y_up;
err = err(:);
err2 = dot(err,err);

MaxIter = 10;
for iter = 1:MaxIter,
%     fprintf(1,'%d...',iter);
    if updateJ,
        om = param(1:3);
        T = param(4:6);
        [~,dYdom,dYdT] = rigid_refmotion(X,om,T,hand);
        JJ = [dYdom, dYdT];
        % compute the approximated Hessian matrix
        JJ2 = JJ'*JJ;
        JJTerr = JJ'*err;
    end;
    H_op = JJ2 + diag(lamda*diag(JJ2));  % JJ2 + lamda*eye(6);
    param_innov = H_op\JJTerr;
    param_up = param + param_innov;
    om_up = param_up(1:3);
    T_up = param_up(4:6);
    Y_up = rigid_refmotion(X,om_up,T_up,hand);
    err = Y - Y_up;
    err = err(:);
    err2_up = dot(err,err);
    if err2_up < err2,
        lamda = lamda/10;
        param = param_up;
        err2 = err2_up;
        updateJ=1;
    else
        lamda = lamda*10;
        updateJ=0;
    end;
    if norm(JJTerr) < bigeps,     % grad = JJTerr;
        break;
    end;
end;

om = param(1:3);
T = param(4:6);
return;



%% Test
% look out the rigid body have symmetric axis

n1 = 10;
X1 = randn(3,n1)*10;
om1 = randn(3,1);
T1 = randn(3,1)*1000;
if n1==3,
    h1 = 1;
else
    h1 = sign(randn);
end;
delta = randn(3,n1);
Y1 = rodrigues(om1)*diag([1,1,h1])*X1+T1 + delta;
[om2,T2,h2] = compute_rigid_refmotion(Y1,X1);
[om3,T3,h3] = compute_rigid_refmotion(Y1,X1,0);
e_op = [om2-om1; T2-T1; h2-h1];
e_nop = [om3-om1; T3-T1; h3-h1];
gain = norm(e_nop)/norm(e_op)
d_init = delta(:)'*delta(:)
d_op = Y1-rigid_refmotion(X1,om2,T2,h2);
d_op =  d_op(:)'*d_op(:)
d_nop = Y1-rigid_refmotion(X1,om3,T3,h3);
d_nop = d_nop(:)'*d_nop(:)

