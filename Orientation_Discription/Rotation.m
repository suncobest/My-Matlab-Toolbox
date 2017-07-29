%% ����ת�������������ת��ϵıȽ�
% ����ת���������ö����Ǹ��壬����ʸ��
% ����ŷ��ת���������������������λ֮�����ͨ����ĳ�ռ����һ��ת����ʵ�֡�
% ���ת��omegaΪ��λʸ��,ת��Ϊfi��
clear, clc
syms w x y z fi ax ay az
% �ڶ�ά�ռ��У���ת�õ�һ�Ƕ�fi���ɱ�ʾ�������þ����ʾ�����ԣ���expm()����ʾ����ָ��,��mpower
% R=[cos(fi),-sin(fi); sin(fi),cos(fi)]

R = simplify( expm(fi*[0,-1;1,0]) )



% ͳһ��ʾ�Ƹ���Ļ�����ת����
% Rx=[1,0,0; 0,cos(ax),-sin(ax); 0,sin(ax),cos(ax)];
% Ry=[cos(ay),0,sin(ay); 0,1,0; -sin(ay),0,cos(ay)];
% Rz=[cos(az),-sin(az),0; sin(az),cos(az),0; 0,0,1];


%  �������ⷽ��r����omega=r/norm(r)��ʸ��r�ĳ���norm(r)=fi
% ����r,�䷴�Գ�����Ϊantir����λʸ��omega�ķ��Գ�����Ϊantiomega��
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
% ��ָ����ʽ�ȼۣ���ת����ĵ�����������������ʽ�������1��������ʽΪ-1��������ת���󣬱����ת�ӷ��ݡ�
Rotmatrix = cos(fi)*eye(3)+(1-cos(fi))*(omega*omega.')+sin(fi)*antiomega

% Rotmatrix = eye(3)+sin(fi)*antiomega+(1-cos(fi))*antiomega^2  



% ��Ԫ��(quaternion)
% ��Ԫ��q = w+xi+yj+zk������[w, x, y, z]����norm(q)=1����w^2+x^2+y^2+z^2=1
% q=cos(fi/2)+sin(fi/2)*omega����������ʸ������Ԫ���˷�����ʹ��ʸ������omega��תfi
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
% % ŷ����ת��ȡ������ڵ�������תZ-Y'-Z"
% r=[sin(b)*cos(a);sin(b)*sin(a);cos(b)];  % ת����ķ�������
% P=[cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1;];  % ��Z�����pusi
% N=[cos(b),0,sin(b);0,1,0;-sin(b),0,cos(b);];  % ��Y'���¶�theta
% R=[cos(c),-sin(c),0;sin(c),cos(c),0;0,0,1;];  % ��Z"����תfi
% r1=cos(c)*[1,0,0;0,1,0;0,0,1;]+(1-cos(c))*(r*r.')+sin(c)*[0,-cos(b),sin(b)*sin(a);cos(b),0,-sin(b)*cos(a);-sin(b)*sin(a),sin(b)*cos(a),0;];  % ����ת������
% r2=P*N*R*N.'*P.';  % ������ת������ϵ��ƿռ���ת���ı任����
% M=r1-r2;   % ȡ��ֵ��֤�����Ƿ����
% simplify(M) 


%% �÷������ң�alfa��beta��garma����ʾת����ķ�����theta��ʾת���ǵĴ�С����������ת�������ľ�����ʽ

% syms alfa beta garma theta % alfa^2+beta^2+garma^2=1
% a=alfa
% b=beta
% c=garma
% 
% r=[a;b;c] % ת����ķ�������
% R=cos(theta)*[1,0,0;0,1,0;0,0,1;]+(1-cos(theta))*r*r.'+sin(theta)*[0,-c,b;c,0,-a;-b,a,0;]
% simplify(R)