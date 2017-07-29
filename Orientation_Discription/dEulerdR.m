% Jacobian matrix of Euler angle over rotation matrix
% The rotation matrix
R=sym('R%d%d',[3 3]);
R=sym(R,'real');

%% dXYXdR
a2 = acos(R(1,1));
a1 = atan2(R(2,1),-R(3,1));
a3 = atan2(R(1,2),R(1,3));
XYX = [a1;a2;a3];
dXYXdR = simplify(jacobian(XYX, R(:)));

%% dXZXdR
a2 = acos(R(1,1));
a1 = atan2(R(3,1),R(2,1));
a3 = atan2(R(1,3),-R(1,2));
XZX = [a1;a2;a3];
dXZXdR = simplify(jacobian(XZX, R(:)));

%% dYXYdR
a2 = acos(R(2,2));
a1 = atan2(R(1,2),R(3,2));
a3 = atan2(R(2,1),-R(2,3));
YXY = [a1;a2;a3];
dYXYdR = simplify(jacobian(YXY, R(:)));

%% dYZYdR
a2 = acos(R(2,2));
a1 = atan2(R(3,2),-R(1,2));
a3 = atan2(R(2,3),R(2,1));
YZY = [a1;a2;a3];
dYZYdR = simplify(jacobian(YZY, R(:)));

%% dZXZdR
a2 = acos(R(3,3));
a1 = atan2(R(1,3),-R(2,3));
a3 = atan2(R(3,1),R(3,2));
ZXZ = [a1;a2;a3];
dZXZdR = simplify(jacobian(ZXZ, R(:)));

%% dZYZdR
a2 = acos(R(3,3));
a1 = atan2(R(2,3),R(1,3));
a3 = atan2(R(3,2),-R(3,1));
ZYZ = [a1;a2;a3];
dZYZdR = simplify(jacobian(ZYZ, R(:)));

%% dXYZdR
a2 = asin(R(1,3));
a1 = atan2(-R(2,3),R(3,3));
a3 = atan2(-R(1,2),R(1,1));
XYZ = [a1;a2;a3];
dXYZdR = simplify(jacobian(XYZ, R(:)));

%% dXZYdR
a2 = -asin(R(1,2));
a1 = atan2(R(3,2),R(2,2));
a3 = atan2(R(1,3),R(1,1));
XZY = [a1;a2;a3];
dXZYdR = simplify(jacobian(XZY, R(:)));

%% dYXZdR
a2 = -asin(R(2,3));
a1 = atan2(R(1,3),R(3,3));
a3 = atan2(R(2,1),R(2,2));
YXZ = [a1;a2;a3];
dYXZdR = simplify(jacobian(YXZ, R(:)));

%% dYZXdR
a2 = asin(R(2,1));
a1 = atan2(-R(3,1),R(1,1));
a3 = atan2(-R(2,3),R(2,2));
YZX = [a1;a2;a3];
dYZXdR = simplify(jacobian(YZX, R(:)));

%% dZXYdR
a2 = asin(R(3,2));
a1 = atan2(-R(1,2),R(2,2));
a3 = atan2(-R(3,1),R(3,3));
ZXY = [a1;a2;a3];
dZXYdR = simplify(jacobian(ZXY, R(:)));

%% dZYXdR
a2 = -asin(R(3,1));
a1 = atan2(R(2,1),R(1,1));
a3 = atan2(R(3,2),R(3,3));
ZYX = [a1;a2;a3];
dZYXdR = simplify(jacobian(ZYX, R(:)));