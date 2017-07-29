function [Qckk,Tckk,JJ] = compute_extrinsic_lm_refine(x_kk,X_kk,Qc_init,Tc_init,hand,fc,cc,kc,alpha_c,MaxIter,thresh_cond)

% compute_extrinsic_lm_refine
%
% Computes the extrinsic parameters attached to a 3D structure X_kk given its projection
% on the image plane x_kk and the intrinsic camera parameters fc, cc and kc.
% Works with planar and non-planar structures.
%
% INPUT: x_kk: Feature locations on the images
%       X_kk: Corresponding grid coordinates
%       fc: Camera focal length
%       cc: Principal point coordinates
%       kc: Distortion coefficients
%       alpha_c: Skew coefficient
%       MaxIter: Maximum number of iterations
%
% OUTPUT: Qckk: 4D quaternion rotation vector attached to the grid positions in space
%         Tckk: 3D translation vector attached to the grid positions in space
%         JJ: the jacobian matrix [dxdQ dxdT].
%
% Method: Computes the normalized point coordinates, then computes the 3D pose
%
% Important functions called within that program:
% normalize_pixel: Computes the normalize image point coordinates.
% trans_quat_mat: transform rotation matrix to quternion vector, or vice versa.
%  See also compute_extrinsic_init2, compute_extrinsic_init.

if nargin < 11,
    thresh_cond = inf;
    if nargin < 10,
        MaxIter = 20;
        if nargin < 9,
            alpha_c = 0;
            if nargin < 8,
                kc = zeros(5,1);
                if nargin < 7,
                    cc = zeros(2,1);
                    if nargin <6,
                        fc = ones(2,1);
                        if nargin < 5,
                            hand = 1;
                            if nargin < 4,
                                error('Need initial extrinsic value with 2D and 3D points!');
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;
bigeps = eps*1e6;

% Initialization:
Qckk = Qc_init;
Tckk = Tc_init;

% Final optimization (minimize the reprojection error in pixel):
% through LM:

param = [Qckk;Tckk];
lamda = 0.001; % set an initial value of the damping factor for the LM
updateJ = 1;
x = project_points_mirror(X_kk,Qckk,Tckk,hand,fc,cc,kc,alpha_c);
ex = x_kk - x;
ex = ex(:);
ex2 = dot(ex,ex);

for iter = 1:MaxIter,
%     fprintf(1,'%d...',iter);
    if updateJ,
        Qckk = param(1:4);
        Tckk = param(5:7);
        [~,dxdQ,dxdT] = project_points_mirror(X_kk,Qckk,Tckk,hand,fc,cc,kc,alpha_c);
        JJ = [dxdQ, dxdT];
        if cond(JJ) > thresh_cond,    % do not refine
            break;
        end;
        % compute the approximated Hessian matrix
        JJ2 = JJ'*JJ;
        JJTex = JJ'*ex;
    end;
    H_lm = JJ2 + diag(lamda*diag(JJ2));  % JJ2 + lamda*eye(7);
    param_innov = H_lm\JJTex;
    param_up = param + param_innov;
    Qckk_up = param_up(1:4);
    Qckk_up = Qckk_up/norm(Qckk_up);
    param_up(1:4) = Qckk_up;
    Tckk_up = param_up(5:7);
    x = project_points_mirror(X_kk,Qckk_up,Tckk_up,hand,fc,cc,kc,alpha_c);
    ex = x_kk - x;
    ex = ex(:);
    ex2_up = dot(ex,ex);
    if ex2_up < ex2,
        lamda = lamda/10;
        param = param_up;
        ex2 = ex2_up;
        updateJ=1;
    else
        lamda = lamda*10;
        updateJ=0;
    end;
    if norm(JJTex) < bigeps,     % grad = JJTex;
        break;
    end;
end;

Qckk = param(1:4);
Tckk = param(5:7);
return;



%% Test
m = 3;
n = 1000;
X = 100*randn(m,n);
om = randn(3,1);
om = om/norm(om)*pi/2*rand;
q = trans_quat_axis(om);
R = trans_quat_mat(q);
T = [20*randn(2,1);1000*rand+500];
hand = sign(randn(1));
f = 1000*rand(2,1);
c = 1000*randn(2,1);
k = 0.4*randn(5,1);
k(2:5) = k(2:5)/10;
alpha = 0.01*randn(1,1);
[x,dxdq,dxdT,dxdf,dxdc,dxdk,dxdalpha,dxdX] = project_points_mirror(X,q,T,hand,f,c,k,alpha);
dX = randn(m,n)/1000;
dom = randn(3,1)/1000;
dq = trans_quat_axis(om+dom)-q;
dT = norm(T)*randn(3,1)/1000;
df = norm(f)*randn(2,1)/1000;
dc = norm(c)*randn(2,1)/100;
dk = norm(k)*randn(5,1)/1000;
dalpha = norm(k)*randn(1,1)/1000;
x2 = project_points_mirror(X+dX,q+dq,T+dT,hand,f+df,c+dc,k+dk,alpha+dalpha);
x_pred = x +reshape(dxdX*dX(:),2,n)+ reshape(dxdq*dq,2,n)+reshape(dxdT*dT,2,n)+reshape(dxdf*df,2,n)+reshape(dxdc*dc,2,n)+reshape(dxdk*dk,2,n)+reshape(dxdalpha*dalpha,2,n);
gain = norm(x2-x)/norm(x2 - x_pred)
%%
[Qc,Tc] = compute_extrinsic_init_mirror(x,X,f,c,k,alpha,hand);
if sum(Qc.*q)<0,
    Qc = -Qc;
end;
err_init =norm([q;T]-[Qc;Tc])
[Qc1,Tc1,JJ1] = compute_extrinsic_lm_refine(x,X,Qc,Tc,hand,f,c,k,alpha,15,1e6);
err_lm =norm([q;T]-[Qc1;Tc1])

