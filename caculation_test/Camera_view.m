clear;clc;echo off;
mesh(peaks);box;
axis square tight
% axis equal;


% x=randi(300)*randn(1,3),y=rot2axis(axis2rot(x))*norm(x),...
% a=randi(300)*randn(1,2),b=axis2rot(rot2axis(a))
% view(x),[az,el]=view;[az,el]
% pause;
% view(y),[az,el]=view;[az,el]
% pause;
% view(a),a
% pause;
% view(b),b

%%

% az_el=randi(300)*randn(1,2)  
% r1=rot2axis(az_el),r2=rot2axis(az_el(1)+180,180-az_el(2))  % 按照定义，r1和r2是同一个轴
% view(az_el),view,pause   % 以下四种情况视角一致
% view(az_el(1)+180,180-az_el(2)),view,pause
% view(r1),view,pause
% view(r2),view,pause



%% 世界坐标系到摄像机坐标系的变换：Xc=Rm*Xw+T
% 
% 其中Xc为摄像机坐标；Xw为世界坐标；Rm为旋转矩阵rotation matrix，Rm的列向量分别
% 表示世界坐标系基矢量在摄像机坐标系下的投影分量；T为世界坐标系原点在摄像机坐标系中的位置。
% 
% 在摄像机坐标系中：y轴作为景深方向，指向前方；z轴朝上；x轴朝右。从世界坐标系到摄像
% 机坐标系的变换，包括旋转加平移。按照动轴旋转（内禀旋转）定义，旋转包括三个欧拉角：
% 首先是绕z轴旋转角度a1，接着是绕x轴旋转角度a2，再者是绕y轴旋转角度a3。其中a1是方
% 位角（Azimuth），a2是俯仰角（Elevation），a2是滚转角（Roll）。即a1=az；另外由于
% 视线沿y轴，-y方向与世界坐标系水平面夹角决定了el，所以a2=-el；最后摄像机坐标系绕
% 视线yc轴正向旋转时，画面正好逆时针旋转，所以画面倾斜角a3=ro。若摄像机不绕视线旋转，
% 则a3取0。
% 
% 统一表示绕各轴的基本旋转矩阵
% Rx=[1,0,0; 0,cos(ax),-sin(ax); 0,sin(ax),cos(ax)];
% Ry=[cos(ay),0,sin(ay); 0,1,0; -sin(ay),0,cos(ay)];
% Rz=[cos(az),-sin(az),0; sin(az),cos(az),0; 0,0,1];
%
% 按照内禀旋转的定义，摄像机坐标系的欧拉旋转为ZXY,即az=a1,ax=a2,ay=a3,Rm=ZXY。
 
% syms a1 a2 a3  % 生成三个欧拉角变量

syms az el ro  

az = randi(300)*randn(1);
el = randi(300)*randn(1);
ro = 0;

a1 = az;
a2 = -el;
a3 = ro;

% % Rm=ZXY
Rm(1, :) = [cos(a1)*cos(a3) - sin(a1)*sin(a2)*sin(a3), -cos(a2)*sin(a1), cos(a1)*sin(a3) + cos(a3)*sin(a1)*sin(a2)];
Rm(2, :) = [cos(a3)*sin(a1) + cos(a1)*sin(a2)*sin(a3), cos(a1)*cos(a2), sin(a1)*sin(a3) - cos(a1)*cos(a3)*sin(a2)];
Rm(3, :) = [-cos(a2)*sin(a3), sin(a2), cos(a2)*cos(a3)];

% Rm = subs( Rm, {a1,a2,a3}, {sym('az'),-sym('el'),sym('ro')} );   % 替换符号变量

T = zeros(3,1);
Xw = [0;-1;0];
Xc = Rm*Xw+T;  



% 在世界坐标系下，Rm作为旋转算子，作用于矢量Xw在得出旋转之后的新位矢。
% 新位矢在世界坐标系下的分量，等于矢量Xw在摄像机坐标系下的投影分量Xc。

view(az,el);
Rt = view 
Rm 


% Xw = Rm.'*(Xc-T); 

K = [fx,0,0;s,fy,0;cx,cy,1]

s*[u,v,1]=[xw,yw,zw,1]*[R';T']*K


% The coordinates [cx cy]
% represent the optical center (the principal point), in pixels. When
% the x and y axis are exactly
% perpendicular, the skew parameter, s, equals 0.
% fx = F*sx
% fy = F*sy
% where 
% F, is the focal length in world units,typically expressed in millimeters.
% [sx, sy] are the number of pixels per world unit in the x and y direction respectively. 
% fx and fy are expressed in pixels.

sx=1/du
sy=1/dv

zc*[u;v;1]=K*[R,T]*[xw;yw;zw;1]

% az=0;
% el=0;
% % Rx*Rz*[0;-1;0]
% R=simplify(R),R=subs(R,[a1,a2,a3],[az,el,0])
% % subs(R,{az,el},{0,pi/2})

