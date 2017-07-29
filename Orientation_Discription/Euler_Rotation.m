%% 有限转动张量与基本旋转组合的比较
% 有限转动张量作用对象是刚体，例如矢量；与坐标系变换区分开，相当于坐标系变换的逆变换

% 统一表示绕各轴的基本旋转矩阵
%     Rx=[1,0,0; 0,cos(ax),-sin(ax); 0,sin(ax),cos(ax)];
%     Ry=[cos(ay),0,sin(ay); 0,1,0; -sin(ay),0,cos(ay)];
%     Rz=[cos(az),-sin(az),0; sin(az),cos(az),0; 0,0,1];


% 用相互独立的进动角(经度)和章动角(维度)表示转动轴的方向，fi表示转动角
syms pusi theta fi  % pusi是进动角(Precession)，theta是章动角(Nutation)，fi是自转角(Rotation)
a=pusi;
b=theta;
c=fi;

% 欧拉旋转采取由外而内的内禀旋转Z-Y'-Z"
r=[sin(b)*cos(a);sin(b)*sin(a);cos(b)];  % 转动轴的方向余弦
P=[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1;];  % 绕Z轴进动pusi
N=[cos(b),0,sin(b);0,1,0;-sin(b),0,cos(b);];  % 绕Y'轴章动theta
R=[cos(c),-sin(c),0;sin(c),cos(c),0;0,0,1;];  % 绕Z"轴自转fi
r1=cos(c)*[1,0,0;0,1,0;0,0,1;]+(1-cos(c))*(r*r.')+sin(c)*[0,-cos(b),sin(b)*sin(a);cos(b),0,-sin(b)*cos(a);-sin(b)*sin(a),sin(b)*cos(a),0;];  % 有限转动张量
r2=P*N*R*N.'*P.';  % 基本旋转矩阵组合的绕空间轴转动的变换矩阵
M=r1-r2;   % 取差值验证两者是否相等
simplify(M) 


%% 用方向余弦（alfa，beta，garma）表示转动轴的方向，用theta表示转动角的大小，给出有限转动张量的矩阵表达式

% syms alfa beta garma theta % alfa^2+beta^2+garma^2=1
% a=alfa
% b=beta
% c=garma
% 
% r=[a;b;c] % 转动轴的方向余弦
% R=cos(theta)*[1,0,0;0,1,0;0,0,1;]+(1-cos(theta))*r*r.'+sin(theta)*[0,-c,b;c,0,-a;-b,a,0;]
% simplify(R)