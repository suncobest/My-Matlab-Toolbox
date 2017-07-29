% transform change rate of Euler angle to angular velocity
% three Euler angle: a1, a2, a3
% change rate of three angle: a1t, a2t, a3t
% 1 denotes intrinsic notation; 0 denotes extrinsic notation;

syms a1 a2 a3 a1t a2t a3t om1 om2 om3 real

rotx=@(x) [1,0,0; 0,cos(x),-sin(x); 0,sin(x),cos(x)];
roty=@(x) [cos(x),0,sin(x); 0,1,0; -sin(x),0,cos(x)];
rotz=@(x) [cos(x),-sin(x),0; sin(x),cos(x),0; 0,0,1];
at = [a1t; a2t; a3t];
aa = [a1; a2; a3];
om = [om1; om2; om3];

%% XYX
Rxyx = simplify(rotx(a1)*roty(a2)*rotx(a3));
OMxyx1 = a1t*Rxyx(1,:)'+a2t*rotx(a3)'*[0;1;0]+a3t*[1;0;0];
OMxyx0 = a1t*[1;0;0]+a2t*rotx(a1)*[0;1;0]+a3t*Rxyx(:,1);
dOMxyx1dat = simplify(jacobian(OMxyx1,at));
dOMxyx1daa = simplify(jacobian(OMxyx1,aa));
dOMxyx0dat = simplify(jacobian(OMxyx0,at));
dOMxyx0daa = simplify(jacobian(OMxyx0,aa));
simplify(OMxyx0-Rxyx*OMxyx1)

datdOMxyx1 = simplify(dOMxyx1dat\eye(3));
datdOMxyx0 = simplify(dOMxyx0dat\eye(3));
datdaa_xyx1 =simplify(jacobian(dOMxyx1dat\om,aa));
datdaa_xyx0 =simplify(jacobian(dOMxyx0dat\om,aa));

%% XZX
Rxzx = simplify(rotx(a1)*rotz(a2)*rotx(a3));
OMxzx1 = a1t*Rxzx(1,:)'+a2t*rotx(a3)'*[0;0;1]+a3t*[1;0;0];
OMxzx0 = a1t*[1;0;0]+a2t*rotx(a1)*[0;0;1]+a3t*Rxzx(:,1);
dOMxzx1dat = simplify(jacobian(OMxzx1,at));
dOMxzx1daa = simplify(jacobian(OMxzx1,aa));
dOMxzx0dat = simplify(jacobian(OMxzx0,at));
dOMxzx0daa = simplify(jacobian(OMxzx0,aa));
simplify(OMxzx0-Rxzx*OMxzx1)

datdOMxzx1 = simplify(dOMxzx1dat\eye(3));
datdOMxzx0 = simplify(dOMxzx0dat\eye(3));
datdaa_xzx1 =simplify(jacobian(dOMxzx1dat\om,aa));
datdaa_xzx0 =simplify(jacobian(dOMxzx0dat\om,aa));

%% YXY
Ryxy = simplify(roty(a1)*rotx(a2)*roty(a3));
OMyxy1 = a1t*Ryxy(2,:)'+a2t*roty(a3)'*[1;0;0]+a3t*[0;1;0];
OMyxy0 = a1t*[0;1;0]+a2t*roty(a1)*[1;0;0]+a3t*Ryxy(:,2);
dOMyxy1dat = simplify(jacobian(OMyxy1,at));
dOMyxy1daa = simplify(jacobian(OMyxy1,aa));
dOMyxy0dat = simplify(jacobian(OMyxy0,at));
dOMyxy0daa = simplify(jacobian(OMyxy0,aa));
simplify(OMyxy0-Ryxy*OMyxy1)

datdOMyxy1 = simplify(dOMyxy1dat\eye(3));
datdOMyxy0 = simplify(dOMyxy0dat\eye(3));
datdaa_yxy1 =simplify(jacobian(dOMyxy1dat\om,aa));
datdaa_yxy0 =simplify(jacobian(dOMyxy0dat\om,aa));

%% YZY
Ryzy = simplify(roty(a1)*rotz(a2)*roty(a3));
OMyzy1 = a1t*Ryzy(2,:)'+a2t*roty(a3)'*[0;0;1]+a3t*[0;1;0];
OMyzy0 = a1t*[0;1;0]+a2t*roty(a1)*[0;0;1]+a3t*Ryzy(:,2);
dOMyzy1dat = simplify(jacobian(OMyzy1,at));
dOMyzy1daa = simplify(jacobian(OMyzy1,aa));
dOMyzy0dat = simplify(jacobian(OMyzy0,at));
dOMyzy0daa = simplify(jacobian(OMyzy0,aa));
simplify(OMyzy0-Ryzy*OMyzy1)

datdOMyzy1 = simplify(dOMyzy1dat\eye(3));
datdOMyzy0 = simplify(dOMyzy0dat\eye(3));
datdaa_yzy1 =simplify(jacobian(dOMyzy1dat\om,aa));
datdaa_yzy0 =simplify(jacobian(dOMyzy0dat\om,aa));

%% ZXZ
Rzxz = simplify(rotz(a1)*rotx(a2)*rotz(a3));
OMzxz1 = a1t*Rzxz(3,:)'+a2t*rotz(a3)'*[1;0;0]+a3t*[0;0;1];
OMzxz0 = a1t*[0;0;1]+a2t*rotz(a1)*[1;0;0]+a3t*Rzxz(:,3);
dOMzxz1dat = simplify(jacobian(OMzxz1,at));
dOMzxz1daa = simplify(jacobian(OMzxz1,aa));
dOMzxz0dat = simplify(jacobian(OMzxz0,at));
dOMzxz0daa = simplify(jacobian(OMzxz0,aa));
simplify(OMzxz0-Rzxz*OMzxz1)

datdOMzxz1 = simplify(dOMzxz1dat\eye(3));
datdOMzxz0 = simplify(dOMzxz0dat\eye(3));
datdaa_zxz1 =simplify(jacobian(dOMzxz1dat\om,aa));
datdaa_zxz0 =simplify(jacobian(dOMzxz0dat\om,aa));

%% ZYZ
Rzyz = simplify(rotz(a1)*roty(a2)*rotz(a3));
OMzyz1 = a1t*Rzyz(3,:)'+a2t*rotz(a3)'*[0;1;0]+a3t*[0;0;1];
OMzyz0 = a1t*[0;0;1]+a2t*rotz(a1)*[0;1;0]+a3t*Rzyz(:,3);
dOMzyz1dat = simplify(jacobian(OMzyz1,at));
dOMzyz1daa = simplify(jacobian(OMzyz1,aa));
dOMzyz0dat = simplify(jacobian(OMzyz0,at));
dOMzyz0daa = simplify(jacobian(OMzyz0,aa));
simplify(OMzyz0-Rzyz*OMzyz1)

datdOMzyz1 = simplify(dOMzyz1dat\eye(3));
datdOMzyz0 = simplify(dOMzyz0dat\eye(3));
datdaa_zyz1 =simplify(jacobian(dOMzyz1dat\om,aa));
datdaa_zyz0 =simplify(jacobian(dOMzyz0dat\om,aa));

%% XYZ
Rxyz = simplify(rotx(a1)*roty(a2)*rotz(a3));
OMxyz1 = a1t*Rxyz(1,:)'+a2t*rotz(a3)'*[0;1;0]+a3t*[0;0;1];
OMxyz0 = a1t*[1;0;0]+a2t*rotx(a1)*[0;1;0]+a3t*Rxyz(:,3);
dOMxyz1dat = simplify(jacobian(OMxyz1,at));
dOMxyz1daa = simplify(jacobian(OMxyz1,aa));
dOMxyz0dat = simplify(jacobian(OMxyz0,at));
dOMxyz0daa = simplify(jacobian(OMxyz0,aa));
simplify(OMxyz0-Rxyz*OMxyz1)

datdOMxyz1 = simplify(dOMxyz1dat\eye(3));
datdOMxyz0 = simplify(dOMxyz0dat\eye(3));
datdaa_xyz1 =simplify(jacobian(dOMxyz1dat\om,aa));
datdaa_xyz0 =simplify(jacobian(dOMxyz0dat\om,aa));

%% XZY
Rxzy = simplify(rotx(a1)*rotz(a2)*roty(a3));
OMxzy1 = a1t*Rxzy(1,:)'+a2t*roty(a3)'*[0;0;1]+a3t*[0;1;0];
OMxzy0 = a1t*[1;0;0]+a2t*rotx(a1)*[0;0;1]+a3t*Rxzy(:,2);
dOMxzy1dat = simplify(jacobian(OMxzy1,at));
dOMxzy1daa = simplify(jacobian(OMxzy1,aa));
dOMxzy0dat = simplify(jacobian(OMxzy0,at));
dOMxzy0daa = simplify(jacobian(OMxzy0,aa));
simplify(OMxzy0-Rxzy*OMxzy1)

datdOMxzy1 = simplify(dOMxzy1dat\eye(3));
datdOMxzy0 = simplify(dOMxzy0dat\eye(3));
datdaa_xzy1 =simplify(jacobian(dOMxzy1dat\om,aa));
datdaa_xzy0 =simplify(jacobian(dOMxzy0dat\om,aa));

%% YXZ
Ryxz = simplify(roty(a1)*rotx(a2)*rotz(a3));
OMyxz1 = a1t*Ryxz(2,:)'+a2t*rotz(a3)'*[1;0;0]+a3t*[0;0;1];
OMyxz0 = a1t*[0;1;0]+a2t*roty(a1)*[1;0;0]+a3t*Ryxz(:,3);
dOMyxz1dat = simplify(jacobian(OMyxz1,at));
dOMyxz1daa = simplify(jacobian(OMyxz1,aa));
dOMyxz0dat = simplify(jacobian(OMyxz0,at));
dOMyxz0daa = simplify(jacobian(OMyxz0,aa));
simplify(OMyxz0-Ryxz*OMyxz1)

datdOMyxz1 = simplify(dOMyxz1dat\eye(3));
datdOMyxz0 = simplify(dOMyxz0dat\eye(3));
datdaa_yxz1 =simplify(jacobian(dOMyxz1dat\om,aa));
datdaa_yxz0 =simplify(jacobian(dOMyxz0dat\om,aa));

%% YZX
Ryzx = simplify(roty(a1)*rotz(a2)*rotx(a3));
OMyzx1 = a1t*Ryzx(2,:)'+a2t*rotx(a3)'*[0;0;1]+a3t*[1;0;0];
OMyzx0 = a1t*[0;1;0]+a2t*roty(a1)*[0;0;1]+a3t*Ryzx(:,1);
dOMyzx1dat = simplify(jacobian(OMyzx1,at));
dOMyzx1daa = simplify(jacobian(OMyzx1,aa));
dOMyzx0dat = simplify(jacobian(OMyzx0,at));
dOMyzx0daa = simplify(jacobian(OMyzx0,aa));
simplify(OMyzx0-Ryzx*OMyzx1)

datdOMyzx1 = simplify(dOMyzx1dat\eye(3));
datdOMyzx0 = simplify(dOMyzx0dat\eye(3));
datdaa_yzx1 =simplify(jacobian(dOMyzx1dat\om,aa));
datdaa_yzx0 =simplify(jacobian(dOMyzx0dat\om,aa));

%% ZXY
Rzxy = simplify(rotz(a1)*rotx(a2)*roty(a3));
OMzxy1 = a1t*Rzxy(3,:)'+a2t*roty(a3)'*[1;0;0]+a3t*[0;1;0];
OMzxy0 = a1t*[0;0;1]+a2t*rotz(a1)*[1;0;0]+a3t*Rzxy(:,2);
dOMzxy1dat = simplify(jacobian(OMzxy1,at));
dOMzxy1daa = simplify(jacobian(OMzxy1,aa));
dOMzxy0dat = simplify(jacobian(OMzxy0,at));
dOMzxy0daa = simplify(jacobian(OMzxy0,aa));
simplify(OMzxy0-Rzxy*OMzxy1)

datdOMzxy1 = simplify(dOMzxy1dat\eye(3));
datdOMzxy0 = simplify(dOMzxy0dat\eye(3));
datdaa_zxy1 =simplify(jacobian(dOMzxy1dat\om,aa));
datdaa_zxy0 =simplify(jacobian(dOMzxy0dat\om,aa));

%% ZYX
Rzyx = simplify(rotz(a1)*roty(a2)*rotx(a3));
OMzyx1 = a1t*Rzyx(3,:)'+a2t*rotx(a3)'*[0;1;0]+a3t*[1;0;0];
OMzyx0 = a1t*[0;0;1]+a2t*rotz(a1)*[0;1;0]+a3t*Rzyx(:,1);
dOMzyx1dat = simplify(jacobian(OMzyx1,at));
dOMzyx1daa = simplify(jacobian(OMzyx1,aa));
dOMzyx0dat = simplify(jacobian(OMzyx0,at));
dOMzyx0daa = simplify(jacobian(OMzyx0,aa));
simplify(OMzyx0-Rzyx*OMzyx1)

datdOMzyx1 = simplify(dOMzyx1dat\eye(3));
datdOMzyx0 = simplify(dOMzyx0dat\eye(3));
datdaa_zyx1 =simplify(jacobian(dOMzyx1dat\om,aa));
datdaa_zyx0 =simplify(jacobian(dOMzyx0dat\om,aa));
