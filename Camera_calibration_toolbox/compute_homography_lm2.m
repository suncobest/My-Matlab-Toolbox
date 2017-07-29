function H = compute_homography_lm2(m,M)

%compute_homography_lm2
%
% H = compute_homography_lm2(m,M)
%
%Computes the planar homography between the point coordinates on the plane (M) and the image
%point coordinates (m). lm is the abbreviation for Leverberg-Marquardt algorithm.
%
%INPUT: m: homogeneous coordinates in the image plane (3xN matrix)
%       M: homogeneous coordinates in the plane in 3D (3xN matrix)
%
%OUTPUT: H: Homography matrix (3x3 homogeneous matrix)
%
%Definition: m ~ H*M where "~" means equal up to a non zero scalar factor.
%
%Method: First computes an initial guess for the homography through quasi-linear method.
%        Then, if the total number of points is larger than 4, optimize the solution by minimizing
%        the reprojection error using Leverberg-Marquardt iteration algorithm(in the least squares sense).



Np = size(m,2);

if Np<4,
    error('It takes at least 4 points to compute homography!');
end;
bigeps = eps*1e6;

if size(m,1)<3,
    m = [m;ones(1,Np)];
else
    m = m ./ (ones(3,1)*m(3,:));
end;

if size(M,1)<3,
    M = [M;ones(1,Np)];
else
    M = M ./ (ones(3,1)*M(3,:));
end;

% Prenormalization of point coordinates (very important): 
% (Similarity)

mc = mean(m,2);   % m的几何中心
mscale = sqrt(2)/mean(sqrt((m(1,:)-mc(1)).^2 + (m(2,:)-mc(2)).^2));
Tm = [mscale 0 -mscale*mc(1); 0 mscale -mscale*mc(2); 0 0 1];
inv_Tm = [1/mscale 0 mc(1) ; 0 1/mscale mc(2); 0 0 1];

Mc = mean(M,2);   % M的几何中心
Mscale = sqrt(2)/mean(sqrt((M(1,:)-Mc(1)).^2 + (M(2,:)-Mc(2)).^2));
TM = [Mscale 0 -Mscale*Mc(1); 0 Mscale -Mscale*Mc(2); 0 0 1];

mn = Tm*m; 
Mn = TM*M;

% mn = Hrem * Mn
% Final homography: H = inv_Tm*Hrem*TM;

% 由于L0*h=0，分解测量矩阵L0=A*B，构造大矩阵L=[B,-I;0,A]，则有L*[h;B*h]=0
% Build the matrix: L=[B(6n,9),-eye(6n);zers(2n,9),A(2n,6n)]

L = sparse(8*Np,6*Np+9);                  % 使用稀疏矩阵，初始化为0
L(1:6*Np,10:end) = -speye(6*Np);          % 填入单位阵I

% 填入B
L(1:6*Np,1:6) = repmat(speye(6),Np,1); 
L(1:6*Np,7:9) = -spdiags(reshape(repmat(reshape(mn(1:2,:), 1,[]), 3,1), [],1), 0,6*Np,6*Np) * repmat(speye(3),2*Np,1);

% 填入A
A = sparse(6*Np+3,2*Np-1);
A0 = reshape(repmat(Mn,2,1), 3,[]);
A(1:3,:) = A0(:,1:end-1);
A = reshape([A(:); Mn(:,end)], 6*Np,2*Np)';
L(6*Np+1:end,10:end) = A;

% L*hhv=0
% To find all the singular values of such a matrix, SVD(FULL(A)) will
% usually perform better than svds(A,MIN(SIZE(A))). 
L= L'*L;
[~,~,V] = svd(full(L));

hh = V(1:9,end);  % 最小奇异值对应的右奇异矢量，且norm(hh)=1
Hrem = reshape(hh,3,3)';   % Hrem为Mn到mn的单应性矩阵

flag = 0;   % check H(9)
if abs(Hrem(9))>bigeps,
    Hrem = Hrem/Hrem(9);
    flag = 1;
end;

% keyboard;

%%% Homography refinement if there are more than 4 points:
if Np > 4,
    % Final refinement:
    hhv = reshape(Hrem',9,1);
    MaxIter = 20;
    lamda = 0.001; % set an initial value of the damping factor for the LM
    updateJ = 1;
    mrep = Hrem * Mn;   % reproject
    mrepn = mrep ./ (ones(3,1)*mrep(3,:));  % 归一化，得到mrepn
    err = mn(1:2,:) - mrepn(1:2,:);
    err = err(:);
    err2 = dot(err,err);
    if flag,
        hhv = hhv(1:8);
        for iter=1:MaxIter,
            % jocobi matrix :  J=jacobian(m(1:2)-mrepn(1:2),hhv);
            % J(1, 4:6) = 0;
            % J(2, 1:3) = 0;
            % J(1, 1:3) = -M'/mrep(3);
            % J(2, 4:6) = -M'/mrep(3);
            % J(1, 7:9) = M'*mrep(1)/mrep(3)^2;
            % J(2, 7:9) = M'*mrep(2)/mrep(3)^2;
            % 经测试，9个参数进行迭代时，海森矩阵退化严重，cond(J2）--->∞，因此选取8个参数
            if updateJ,
                J = zeros(2*Np,8);
                MMM = (Mn ./ (ones(3,1)*mrep(3,:)));
                MMM1 = (ones(3,1)*mrepn(1,:)) .* MMM;
                MMM2 = (ones(3,1)*mrepn(2,:)) .* MMM;
                J(1:2:2*Np,1:3) = -MMM';
                J(2:2:2*Np,4:6) = -MMM';
                J(1:2:2*Np,7:8) = MMM1(1:2,:)';
                J(2:2:2*Np,7:8) = MMM2(1:2,:)';
                % compute the approximated Hessian matrix J2
                J2 = J'*J;
                JTerr = J'*err;
            end;
            J2_lm = J2 + lamda*eye(8);     %or  J2_lm = J2 + diag(lamda*diag(J2));
            dhhv = -J2_lm\JTerr;
            hhv_up = hhv + dhhv;
            Hrem = reshape([hhv_up;1],3,3)';
            % reproject
            mrep = Hrem * Mn;
            mrepn = mrep ./ (ones(3,1)*mrep(3,:));  % 归一化，得到mrepn
            err = mn(1:2,:) - mrepn(1:2,:);
            err = err(:);
            err2_up = dot(err,err);
            if err2_up < err2,
                lamda = lamda/10;
                hhv = hhv_up;
                err2 = err2_up;
                updateJ=1;
            else
                lamda = lamda*10;
                updateJ=0;
            end;
%             fprintf(1,'%d...',iter);
            if norm(JTerr) < bigeps,      % grad = JTerr
                break;
            end;
        end;
        Hrem = reshape([hhv;1],3,3)';
    else    
        for iter=1:MaxIter,
            if updateJ,
                J = zeros(2*Np,9);
                MMM = (Mn ./ (ones(3,1)*mrep(3,:)));
                MMM1 = (ones(3,1)*mrepn(1,:)) .* MMM;
                MMM2 = (ones(3,1)*mrepn(2,:)) .* MMM;
                J(1:2:2*Np,1:3) = -MMM';
                J(2:2:2*Np,4:6) = -MMM';
                J(1:2:2*Np,7:9) = MMM1';
                J(2:2:2*Np,7:9) = MMM2';
                % compute the approximated Hessian matrix J2
                J2 = J'*J;
                JTerr = J'*err;
            end;
            J2_lm = J2 + lamda*eye(9);     %or  J2_lm = J2 + diag(lamda*diag(J2));
            dhhv = -J2_lm\JTerr;
            hhv_up = hhv + dhhv;
            Hrem = reshape(hhv_up,3,3)';
            % reproject
            mrep = Hrem * Mn;
            mrepn = mrep ./ (ones(3,1)*mrep(3,:));  % 归一化，得到mrepn
            err = mn(1:2,:) - mrepn(1:2,:);
            err = err(:);
            err2_up = dot(err,err);
            if err2_up < err2,
                lamda = lamda/10;
                hhv = hhv_up;
                err2 = err2_up;
                updateJ=1;
            else
                lamda = lamda*10;
                updateJ=0;
            end;
%             fprintf(1,'%d...',iter);
            if norm(JTerr) < bigeps,      % grad = JTerr
                break;
            end; 
        end;
        Hrem = reshape(hhv,3,3)';
    end;
end;

% Final homography:   mn = Hrem * Mn
H = inv_Tm*Hrem*TM;
if abs(H(9))>bigeps,
    H = H/H(9);
end;
return;



%% test
n=10;
H=randn(3);
H=H/H(9);
X=10*randn(2,n);
x=H*[X;ones(1,n)];
x=x(1:2,:)./(ones(2,1)*x(3,:));
x1=x+randn(2,n)/100;
Hest=compute_homography_lm2(x1,X);
norm(H(:)-Hest(:))
