%% 有限转动张量与基本旋转组合的比较
% 有限转动张量作用对象是刚体，例如矢量
% 根据欧拉转动定理，刚体的任意两个方位之间可以通过绕某空间轴的一次转动来实现。
% 设此转轴omega为单位矢量,转角为fi。
clear, clc
syms w x y z fi ax ay az
% 在二维空间中，旋转用单一角度fi即可表示。除了用矩阵表示还可以，用expm()来表示矩阵指数,见mpower
% R=[cos(fi),-sin(fi); sin(fi),cos(fi)]

R = simplify( expm(fi*[0,-1;1,0]) )



% 统一表示绕各轴的基本旋转矩阵
% Rx=[1,0,0; 0,cos(ax),-sin(ax); 0,sin(ax),cos(ax)];
% Ry=[cos(ay),0,sin(ay); 0,1,0; -sin(ay),0,cos(ay)];
% Rz=[cos(az),-sin(az),0; sin(az),cos(az),0; 0,0,1];


%  给定任意方向r，则omega=r/norm(r)。矢量r的长度norm(r)=fi
% 对于r,其反对称张量为antir；单位矢量omega的反对称张量为antiomega；
% r = 100*rand(3,1);
r = [x; y; z];
fi = norm(r);
omega = r/fi;
antir = [0, -r(3), r(2); r(3), 0, -r(1); -r(2), r(1), 0];
antiomega = [0, -omega(3), omega(2); omega(3), 0, -omega(1); -omega(2), omega(1), 0];
Rot = simplify( expm( fi*antiomega ) )
% Rot = simplify(expm(antir))
Rot = expm(antir)

Rx = simplify( expm(ax*[0,0,0; 0,0,-1; 0,1,0]) )  % r=[1 0 0]
Ry = simplify( expm(ay*[0,0,1; 0,0,0; -1,0,0]) )  % r=[0 1 0]
Rz = simplify( expm(az*[0,-1,0; 1,0,0; 0,0,0]) )  % r=[0 0 1]

% Rodrigues' rotation formula 
% 与指数形式等价，旋转矩阵的第三主不变量即行列式必须等于1，若行列式为-1，则不是旋转矩阵，变成旋转加反演。
Rotmatrix = cos(fi)*eye(3)+(1-cos(fi))*(omega*omega.')+sin(fi)*antiomega

% Rotmatrix = eye(3)+sin(fi)*antiomega+(1-cos(fi))*antiomega^2  



% 四元数(quaternion)
% 四元数q = w+xi+yj+zk，其中[w, x, y, z]满足norm(q)=1。即w^2+x^2+y^2+z^2=1
% q=cos(fi/2)+sin(fi/2)*omega作用于任意矢量（四元数乘法），使得矢量绕轴omega旋转fi
% % Rot = Rotmatrix = Rmat

% quat2mat
Rmat(1,:) = [1-2*y^2-2*z^2, 2*(x*y-w*z), 2*(x*z+w*y), 0];
Rmat(2,:) = [2*(x*y+w*z), 1-2*x^2-2*z^2, 2*(y*z-w*x), 0];
Rmat(3,:) = [2*(x*z-w*y), 2*(w*x+y*z), 1-2*x^2-2*y^2, 0];
Rmat(4,:) = [0, 0, 0, 1]

% mat2quat
w = sqrt(sum( diag(Rmat) ))/2;
x = (Rmat(3,2)- Rmat(2,3))/(4*w);
y = (Rmat(1,3)- Rmat(3,1))/(4*w);
z = (Rmat(2,1)- Rmat(1,2))/(4*w);


% 
% % 欧拉旋转采取由外而内的内禀旋转Z-Y'-Z"
% r=[sin(b)*cos(a);sin(b)*sin(a);cos(b)];  % 转动轴的方向余弦
% P=[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1;];  % 绕Z轴进动pusi
% N=[cos(b),0,sin(b);0,1,0;-sin(b),0,cos(b);];  % 绕Y'轴章动theta
% R=[cos(c),-sin(c),0;sin(c),cos(c),0;0,0,1;];  % 绕Z"轴自转fi
% r1=cos(c)*[1,0,0;0,1,0;0,0,1;]+(1-cos(c))*(r*r.')+sin(c)*[0,-cos(b),sin(b)*sin(a);cos(b),0,-sin(b)*cos(a);-sin(b)*sin(a),sin(b)*cos(a),0;];  % 有限转动张量
% r2=P*N*R*N.'*P.';  % 基本旋转矩阵组合的绕空间轴转动的变换矩阵
% M=r1-r2;   % 取差值验证两者是否相等
% simplify(M) 


%% 用方向余弦（alfa，beta，garma）表示转动轴的方向，用theta表示转动角的大小，给出有限转动张量的矩阵表达式

% syms alfa beta garma theta % alfa^2+beta^2+garma^2=1
% a=alfa
% b=beta
% c=garma
% 
% r=[a;b;c] % 转动轴的方向余弦
% R=cos(theta)*[1,0,0;0,1,0;0,0,1;]+(1-cos(theta))*r*r.'+sin(theta)*[0,-c,b;c,0,-a;-b,a,0;]
% simplify(R)