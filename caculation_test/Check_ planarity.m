%% Check for planarity of the structure:
%keyboard;
% 椭球面
% [X,Y,Z]=ellipsoid(3,2,1,3,2,1,50);
% X_kk=convex_hull([X(:),Y(:),Z(:)]');



% 椭球体
% X_kk=randn(3,10000);  % 球形（立方）正态分布的点阵

X_kk=2*(rand(3,10000)-0.5);    % 立方均匀分布的点阵
vecNorm = @(X) sqrt(sum(X.^2));
X_kk=X_kk(:,vecNorm(X_kk)<=1);   % 球形均匀分布的点阵


Np = size(X_kk,2);

X_kk=diag([3 2 1])*X_kk;   % 拉伸
omega=randn(3,1);
R=rodrigues(omega);
T=10*rand(3,1);
X_kk=R*X_kk+T*ones(1,Np);  % 旋转并平移

X_mean = mean(X_kk')';  % 均值（重心位置）

Y = X_kk - (X_mean*ones(1,Np));  % （误差函数）各点相对重心的位置，X_mean*ones(1,Np)=X_mean(:,ones(1,Np))

YY = Y*Y'; % 协方差矩阵
% 对于三维情形来说，协方差矩阵与转动惯量密切相关。若Np个点的总质量m=1，则转动惯量I=(trace(Y'*Y)*eye(3)-Y*Y')/Np。
% 其中trace(Y'*Y)表示矩阵Y'*Y的迹，等于矩阵Y的所有元素平方和，trace(Y'*Y)=Y(:)'*Y(:)。


[U,S,V] = svd(YY); % 主成分分析PCA
% U和V都是正交矩阵，相当于旋转，V的列向量为旋转后坐标系的基矢量，即主轴的方向矢量。
% 对角矩阵S的对角元素为x^2,y^2和z^2分别对所有点进行求和（积分）后，从大到小排列。
% 相对于主轴的主转动惯量为I=(trace(S)*eye(3)-S)/Np

% YY*V=V*S

% 对于二维情形，相对于主轴的惯性矩为I=(trace(S)*eye(2)-S)/Np
% 对于椭圆x^2/a^2+y^2/b^2<=1来说，其面积为pi*a*b；
% x^2对面积积分得到Iy=pi*a^3*b/4；所以4*Iy/Np=a^2，Iy为椭圆对x轴的惯性矩
% 因此椭圆的半轴长为：
% semi_axes=2*sqrt(diag(S)/Np);


r = S(3,3)/S(2,2);  % 最小奇异值与第二小奇异值的比值

% 画出点阵
figure('Name','ellipsoid','Color','w')
hold on,grid on

plot3(X_kk(1,:),X_kk(2,:),X_kk(3,:),'k.')

axis equal

plot3(X_mean(1),X_mean(2),X_mean(3),'r+')

% 画出三条半主轴

% 对于椭球体x^2/a^2+y^2/b^2+z^2/c^2<=1来说，其体积为4*pi*a*b*c/3；
% x^2对体积积分得到4*pi*a^3*b*c/15；所以5*S(1,1)/Np=a^2。
% 因此椭球的半轴长为：
semi_axes=sqrt(5*diag(S)/Np);

xxx=X_mean+semi_axes(1)*V(:,1);
yyy=X_mean+semi_axes(2)*V(:,2);
zzz=X_mean+semi_axes(3)*V(:,3);

% line([X_mean(1) xxx(1)],[X_mean(2),xxx(2)],[X_mean(3),xxx(3)],'color','r');
% line([X_mean(1) yyy(1)],[X_mean(2),yyy(2)],[X_mean(3),yyy(3)],'color','g');
% line([X_mean(1) zzz(1)],[X_mean(2),zzz(2)],[X_mean(3),zzz(3)],'color','b');

arrow3(X_mean',xxx','r',1,2);
arrow3(X_mean',yyy','g',1,2);
arrow3(X_mean',zzz','b',1,2);

%%
if (r < 1e-3)|(Np < 4), %1e-3, %1e-4, %norm(X_kk(3,:)) < eps, % Test of planarity

   %fprintf(1,'Planar structure detected: r=%f\n',r);

   % Transform the plane to bring it in the Z=0 plane:

   R_transform = V';

   % 酉矩阵V=[v1,v2,v3]，其中v1，v2，v3分别为共面点阵X_kk的三个主方向，且v3为法方向。
   % R_transform=[v1t;v2t;v3t];

   %norm(R_transform(1:2,3))

   if norm(R_transform(1:2,3)) < 1e-6,   % check if R_transform(:,3)==[0;0;1]
      R_transform = eye(3);
   end;

   % 若R_transform(:,3)==[0;0;1]，则v3=[0;0;1]为z方向

   if det(R_transform) < 0, R_transform = -R_transform; end;  % 保证R_transform为旋转，而不是反射

	T_transform = -(R_transform)*X_mean;

	X_new = R_transform*X_kk + T_transform*ones(1,Np);  % X_new为随体坐标（以重心X_mean为原点，主方向v1，v2，v3为基矢量）

    %  X_kk= V*X_new + X_mean*ones(1,Np); （显然X_new(3,:)=0）
