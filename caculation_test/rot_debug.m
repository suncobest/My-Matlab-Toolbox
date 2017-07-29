clear;clc;echo off;
mesh(peaks);
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

% r=randi(300)*randn(1,2),r1=rot2axis(r),r2=rot2axis(r(1)+180,180-r(2))
% view(r1)
% view(r),view
% view(r(1)+180,180-r(2)),view

%% 世界坐标系到摄像机坐标系的变换：Xc=R*Xw+T

% 在摄像机坐标系中：y轴作为景深方向，指向前方；z轴朝上；x轴朝右。从世界坐标系到摄像
% 机坐标系的变换，包括旋转加平移。按照动轴旋转（内禀旋转）定义，旋转包括三个欧拉角：
% 首先是绕z轴旋转角度a1，接着是绕x轴旋转角度a2，再者是绕y轴旋转角度a3。其中a1是方
% 位角（Azimuth），a2是俯仰角（Elevation），a2是滚转角（Roll）。由于视线沿y轴，-y
% 方向与世界坐标系水平面夹角决定了el，所以el=-a2。假设摄像机不绕视线旋转，则a3取0。
% 摄像机坐标系绕视线y轴旋转时，画面正好逆时针旋转，所以画面倾斜角rl=a3。

% 按照内禀旋转的定义，摄像机坐标系的欧拉旋转为ZXY。
% T为世界坐标系原点在摄像机坐标系中的位置。




a='xyz';
j=1;
生成1:3的全排列
syms a1 a2 a3  % 生成三个欧拉角变量
% 统一表示绕各轴的基本旋转矩阵
Rx=[1,0,0; 0,cos(ax),-sin(ax); 0,sin(ax),cos(ax)];
Ry=[cos(ay),0,sin(ay); 0,1,0; -sin(ay),0,cos(ay)];
Rz=[cos(az),-sin(az),0; sin(az),cos(az),0; 0,0,1];

eval(['a',a(i(j,1)),'=a1']);
eval(['a',a(i(j,2)),'=a2']);
R=eval(['R',a(i(j,1)),'*R',a(i(j,2))]);
eval(['a',a(i(j,3)),'=a3']);
R=eval(['R*R',a(i(j,3))]);

disp(j); disp([a(i(j,:)),'=']); disp(R);
j=j+1;
pause
end





% Rx=a2; Ry=a3; Rz=a1;

% Rx=[1,0,0; 0,cos(Rx),-sin(Rx); 0,sin(Rx),cos(Rx)];
% Ry=[cos(Ry),0,sin(Ry); 0,1,0; -sin(Ry),0,cos(Ry)];
% Rz=[cos(Rz),-sin(Rz),0; sin(Rz),cos(Rz),0; 0,0,1];

% R=Rz.'*Rx.'*Ry.';

% Xc=R*Xw+T;
% a1=az;
% a2=-el;

% az=0;
% el=0;
% % Rx*Rz*[0;-1;0]
% R=simplify(R),R=subs(R,[a1,a2,a3],[az,el,0])
% % subs(R,{az,el},{0,pi/2})

