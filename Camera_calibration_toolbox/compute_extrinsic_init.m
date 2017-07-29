function [omckk,Tckk,Rckk] = compute_extrinsic_init(x_kk,X_kk,fc,cc,kc,alpha_c)

%compute_extrinsic
%
%[omckk,Tckk,Rckk] = compute_extrinsic_init(x_kk,X_kk,fc,cc,kc,alpha_c)
%
%Computes the extrinsic parameters attached to a 3D structure X_kk given its projection
%on the image plane x_kk and the intrinsic camera parameters fc, cc and kc.
%Works with planar and non-planar structures.
%
%INPUT: x_kk: Feature locations on the images
%       X_kk: Corresponding grid coordinates
%       fc: Camera focal length
%       cc: Principal point coordinates
%       kc: Distortion coefficients
%       alpha_c: Skew coefficient
%
%OUTPUT: omckk: 3D rotation vector attached to the grid positions in space
%        Tckk: 3D translation vector attached to the grid positions in space
%        Rckk: 3D rotation matrices corresponding to the omc vectors
%
%Method: Computes the normalized point coordinates, then computes the 3D pose
%
%Important functions called within that program:
%
%normalize_pixel: Computes the normalize image point coordinates.
%
%pose3D: Computes the 3D pose of the structure given the normalized image projection.
%
%project_points.m: Computes the 2D image projections of a set of 3D points



if nargin < 6,
    alpha_c = 0;
    if nargin < 5,
        kc = zeros(5,1);
        if nargin < 4,
            cc = zeros(2,1);
            if nargin < 3,
                fc = ones(2,1);
                if nargin < 2,
                    error('Need 2D projections and 3D points (in compute_extrinsic.m)');
                end;
            end;
        end;
    end;
end;

% Compute the normalized coordinates:
xn = normalize_pixel(x_kk,fc,cc,kc,alpha_c);

Np = size(xn,2);

% Check for planarity of the structure:

X_mean = mean(X_kk,2);
Y = X_kk - X_mean*ones(1,Np);
YY = Y*Y';
% Э�������
% ������ά������˵��Э���������ת������������ء���Np�����������m=1��
% ��ת������I=(tr(Y'*Y)*eye(3)-Y*Y')/Np��
% ����tr(Y'*Y)��ʾ����Y'*Y�ļ������ھ���Y������Ԫ��ƽ���ͣ�tr(Y'*Y)=Y(:)'*Y(:)��

% ���ɷַ���PCA�� YY*V=V*S��V��������Ϊ��ϵ�����ᣬS=X_new*X_new',����X_new=V'*(X_kk-X_mean)
[~,S,V] = svd(YY);
% ��С����ֵ��ڶ�С����ֵ�ı�ֵ
r = S(3,3)/S(2,2);

if r<1e-3||Np<4, %norm(X_kk(3,:)) < eps, % Test of planarity
    %fprintf(1,'Planar structure detected: r=%f\n',r);
    % Transform the plane to bring it in the Z=0 plane:
    
    R_trans = V';
    
    % �Ͼ���V=[v1,v2,v3]������v1��v2��v3�ֱ�Ϊ�������X_kk��������������v3Ϊ������
    % ��R_trans(:,3)==[0;0;1]����v3=[0;0;1]Ϊz����
    if norm(R_trans(1:2,3)) < 1e-6,   % check if R_trans(:,3)==[0;0;1]
        R_trans = eye(3);
    end;
    
    % ��֤R_transΪ��ת����������ת+����
    if det(R_trans) < 0,
        R_trans = -R_trans;
    end;
    
    % Compute the planar homography:
    % X_newΪ�������꣨������X_meanΪԭ�㣬������v1��v2��v3Ϊ��ʸ����
    %  X_kk= V*X_new + X_mean*ones(1,Np); ����ȻX_new(3,:)=0��
    X_new = R_trans*Y;    
    H = compute_homography_lm(xn,X_new(1:2,:));  % �����һ��������������ĵ�Ӧ�Ծ���
    
    % H=[h1 h2 h3]=aK[r1 r2 t],aΪ����ֵ����Ϊxn=[Xc/Zc;Yc/Zc;1],����K=eye(3),H=[h1 h2 h3]=[r1 r2 t]/Zc
    % De-embed the motion parameters from the homography:
    sc = mean([norm(H(:,1));norm(H(:,2))]);  % ���h1��h2��ƽ������
    H = H/sc; % ��H���й�һ��
    if H(9)<0,
        H=-H;
    end;
    
    u1 = H(:,1);
    u1 = u1 / norm(u1);
    u2 = H(:,2) - dot(u1,H(:,2)) * u1;
    u2 = u2 / norm(u2);
    u3 = cross(u1,u2);
    Rckk = [u1 u2 u3];
    Tckk = H(:,3);
 
    % Because X_new = R_trans*X_kk + T_trans�� if Xc = Rckk * X_new + Tckk,
    % then Xc = Rckk * R_trans * X_kk + Rckk* T_trans + Tckk
    T_trans = -R_trans*X_mean;
    Tckk = Tckk + Rckk* T_trans;
    Rckk = Rckk * R_trans;    
    omckk = rodrigues(Rckk);
        
else
    %fprintf(1,'Non planar structure detected: r=%f\n',r);
    
    % Computes an initial guess for extrinsic parameters (works for general 3d structure, not planar!!!):
    % The DLT method is applied here!!
    % xn=[r1 r1 r3 t]*[X_kk;1]/Zc
    % ����X_kkΪ���������ά�㣬xnΪ��Zc=1������㣬���߹��ɷǵ�Ӧӳ�䣬
    % ������ֱ�����Ա任DLT��ϵ������[r1 r1 r3 t]/Zc
    
    J = zeros(2*Np,12);
    
    xX = (ones(3,1)*xn(1,:)).*X_kk;
    yX = (ones(3,1)*xn(2,:)).*X_kk;
    
    J(1:2:end,[1 4 7]) = -X_kk';
    J(2:2:end,[2 5 8]) = X_kk';
    J(1:2:end,[3 6 9]) = xX';
    J(2:2:end,[3 6 9]) = -yX';
    J(1:2:end,12) = xn(1,:)';
    J(2:2:end,12) = -xn(2,:)';
    J(1:2:end,10) = -ones(Np,1);
    J(2:2:end,11) = ones(Np,1);
    
    JJ = J'*J;
    [~,~,V] = svd(JJ);
    
    RR = reshape(V(1:9,12),3,3);  % V(1:9,12)=a[r1;r2;r3],aΪ����ʵ��
    
    if det(RR) < 0,
        V(:,12) = -V(:,12);
        RR = -RR;
    end;
    
    [Ur,~,Vr] = svd(RR);
    
    Rckk = Ur*Vr';
    
    sc = norm(V(1:9,12)) / norm(Rckk(:));
    Tckk = V(10:12,12)/sc;
    
    omckk = rodrigues(Rckk);
    
end;

