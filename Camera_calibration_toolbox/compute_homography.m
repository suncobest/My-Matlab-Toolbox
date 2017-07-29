function [H,Hnorm,inv_Hnorm] = compute_homography(m,M);

%compute_homography
%
%[H,Hnorm,inv_Hnorm] = compute_homography(m,M)
%
%Computes the planar homography between the point coordinates on the plane (M) and the image
%point coordinates (m).
%
%INPUT: m: homogeneous coordinates in the image plane (3xN matrix)
%       M: homogeneous coordinates in the plane in 3D (3xN matrix)
%
%OUTPUT: H: Homography matrix (3x3 homogeneous matrix)
%        Hnorm: Normalization matrix used on the points before homography computation
%               (useful for numerical stability is points in pixel coordinates)
%        inv_Hnorm: The inverse of Hnorm
%
%Definition: m ~ H*M where "~" means equal up to a non zero scalar factor.
%
%Method: First computes an initial guess for the homography through quasi-linear method.
%        Then, if the total number of points is larger than 4, optimize the solution by minimizing
%        the reprojection error (in the least squares sense).
%
%
%Important functions called within that program:
%
%comp_distortion_oulu: Undistorts pixel coordinates.
%
%compute_homography.m: Computes the planar homography between points on the grid in 3D, and the image plane.
%
%project_points.m: Computes the 2D image projections of a set of 3D points, and also returns te Jacobian
%                  matrix (derivative with respect to the intrinsic and extrinsic parameters).
%                  This function is called within the minimization loop.




Np = size(m,2);

if size(m,1)<3,
   m = [m;ones(1,Np)];
end;

if size(M,1)<3,
   M = [M;ones(1,Np)];
end;

m = m ./ (ones(3,1)*m(3,:));
M = M ./ (ones(3,1)*M(3,:));

% Prenormalization of point coordinates (very important): 
% (Affine normalization)

ax = m(1,:);
ay = m(2,:);

mxx = mean(ax); % m的几何中心
myy = mean(ay);
ax = ax - mxx;  % m的x误差（相对几何中心的位置）
ay = ay - myy;  % m的y误差

scxx = mean(abs(ax));  % 平均绝对误差
scyy = mean(abs(ay));

Hnorm = [1/scxx 0 -mxx/scxx;0 1/scyy -myy/scyy;0 0 1];
inv_Hnorm = [scxx 0 mxx ; 0 scyy myy; 0 0 1];

% Compute the homography between m and mn:两者之间的单应性矩阵为Hnorm
% Hnorm为点阵m到相对误差点阵mn的变换矩阵（为单应性的）

mn = Hnorm*m;  % 相对误差 =  误差/平均绝对误差：mnx=(mx-mxx)/scxx,mny=(my-myy)/scyy


% Build the matrix:

L = zeros(2*Np,9);

L(1:2:2*Np,1:3) = M';
L(2:2:2*Np,4:6) = M';
L(1:2:2*Np,7:9) = -((ones(3,1)*mn(1,:)).* M)';
L(2:2:2*Np,7:9) = -((ones(3,1)*mn(2,:)).* M)';

% L*hh=0

if Np > 4,
    L = L'*L;
end;

[U,S,V] = svd(L);

hh = V(:,9);  % 最小奇异值对应的右奇异矢量
hh = hh/hh(9); % 规一化

Hrem = reshape(hh,3,3)';   % Hrem为M到mn的单应性矩阵
%Hrem = Hrem / Hrem(3,3);


% Final homography:

H = inv_Hnorm*Hrem;

if 0,
   m2 = H*M;
   m2 = [m2(1,:)./m2(3,:) ; m2(2,:)./m2(3,:)];
   merr = m(1:2,:) - m2;
end;

%keyboard;
 
%%% Homography refinement if there are more than 4 points:

if Np > 4,
   
   % Final refinement:
   hhv = reshape(H',9,1); 
   hhv = hhv(1:8);
   
   % hhv是H的行向量排成一列的前8个元素,H的第9个元素为1
   
   for iter=1:10,
%        fprintf(1,'%d...',iter);
       % reproject
       mrep = H * M;
       
       J = zeros(2*Np,8);
       
       MMM = (M ./ (ones(3,1)*mrep(3,:)));   % 归一化，保证H*MMM的第三行全为1
       
       J(1:2:2*Np,1:3) = -MMM';
       J(2:2:2*Np,4:6) = -MMM';
       
       mrep = mrep ./ (ones(3,1)*mrep(3,:));  % 归一化，保证mrep的第三行全为1
       
       m_err = m(1:2,:) - mrep(1:2,:);
       m_err = m_err(:);
       
       MMM2 = (ones(3,1)*mrep(1,:)) .* MMM;
       MMM3 = (ones(3,1)*mrep(2,:)) .* MMM;
       
       J(1:2:2*Np,7:8) = MMM2(1:2,:)';
       J(2:2:2*Np,7:8) = MMM3(1:2,:)';
       
       % J*hhv = -mrep(1:2,:)(:)
       hh_innov  = inv(J'*J)*J'*m_err;
       
       hhv_up = hhv - hh_innov;
       
       H_up = reshape([hhv_up;1],3,3)';
       
       %norm(m_err)
       %norm(hh_innov)
       
       hhv = hhv_up;
       H = H_up;
%        if norm(J'*m_err) <  eps*1e6,      % grad = JTerr
%            break;
%        end;
   end;
      
end;

return;

%test of Jacobian

mrep = H*M;
mrep = mrep ./ (ones(3,1)*mrep(3,:));

m_err = mrep(1:2,:) - m(1:2,:);
figure(8);
plot(m_err(1,:),m_err(2,:),'r+');
std(m_err')
