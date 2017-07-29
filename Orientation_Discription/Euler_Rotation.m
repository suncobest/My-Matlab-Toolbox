%% ����ת�������������ת��ϵıȽ�
% ����ת���������ö����Ǹ��壬����ʸ����������ϵ�任���ֿ����൱������ϵ�任����任

% ͳһ��ʾ�Ƹ���Ļ�����ת����
%     Rx=[1,0,0; 0,cos(ax),-sin(ax); 0,sin(ax),cos(ax)];
%     Ry=[cos(ay),0,sin(ay); 0,1,0; -sin(ay),0,cos(ay)];
%     Rz=[cos(az),-sin(az),0; sin(az),cos(az),0; 0,0,1];


% ���໥�����Ľ�����(����)���¶���(ά��)��ʾת����ķ���fi��ʾת����
syms pusi theta fi  % pusi�ǽ�����(Precession)��theta���¶���(Nutation)��fi����ת��(Rotation)
a=pusi;
b=theta;
c=fi;

% ŷ����ת��ȡ������ڵ�������תZ-Y'-Z"
r=[sin(b)*cos(a);sin(b)*sin(a);cos(b)];  % ת����ķ�������
P=[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1;];  % ��Z�����pusi
N=[cos(b),0,sin(b);0,1,0;-sin(b),0,cos(b);];  % ��Y'���¶�theta
R=[cos(c),-sin(c),0;sin(c),cos(c),0;0,0,1;];  % ��Z"����תfi
r1=cos(c)*[1,0,0;0,1,0;0,0,1;]+(1-cos(c))*(r*r.')+sin(c)*[0,-cos(b),sin(b)*sin(a);cos(b),0,-sin(b)*cos(a);-sin(b)*sin(a),sin(b)*cos(a),0;];  % ����ת������
r2=P*N*R*N.'*P.';  % ������ת������ϵ��ƿռ���ת���ı任����
M=r1-r2;   % ȡ��ֵ��֤�����Ƿ����
simplify(M) 


%% �÷������ң�alfa��beta��garma����ʾת����ķ�����theta��ʾת���ǵĴ�С����������ת�������ľ�����ʽ

% syms alfa beta garma theta % alfa^2+beta^2+garma^2=1
% a=alfa
% b=beta
% c=garma
% 
% r=[a;b;c] % ת����ķ�������
% R=cos(theta)*[1,0,0;0,1,0;0,0,1;]+(1-cos(theta))*r*r.'+sin(theta)*[0,-c,b;c,0,-a;-b,a,0;]
% simplify(R)