%% 列出12种欧拉角惯例

% clear all
syms a1 a2 a3 real

%% --------------------------------------------------------------------------------------------

XYX(1, :) = [cos(a2), sin(a2)*sin(a3), cos(a3)*sin(a2)];
XYX(2, :) = [sin(a1)*sin(a2), cos(a1)*cos(a3) - cos(a2)*sin(a1)*sin(a3), - cos(a1)*sin(a3) - cos(a2)*cos(a3)*sin(a1)];
XYX(3, :) = [-cos(a1)*sin(a2), cos(a3)*sin(a1) + cos(a1)*cos(a2)*sin(a3), cos(a1)*cos(a2)*cos(a3) - sin(a1)*sin(a3)];
dXYXda = simplify(jacobian(XYX(:), [a1 a2 a3]));

%% --------------------------------------------------------------------------------------------

XZX(1, :) = [cos(a2), -cos(a3)*sin(a2), sin(a2)*sin(a3)];
XZX(2, :) = [cos(a1)*sin(a2), cos(a1)*cos(a2)*cos(a3) - sin(a1)*sin(a3), - cos(a3)*sin(a1) - cos(a1)*cos(a2)*sin(a3)];
XZX(3, :) = [sin(a1)*sin(a2), cos(a1)*sin(a3) + cos(a2)*cos(a3)*sin(a1), cos(a1)*cos(a3) - cos(a2)*sin(a1)*sin(a3)];
dXZXda = simplify(jacobian(XZX(:), [a1 a2 a3]));

%% --------------------------------------------------------------------------------------------

YXY(1, :) = [cos(a1)*cos(a3) - cos(a2)*sin(a1)*sin(a3), sin(a1)*sin(a2), cos(a1)*sin(a3) + cos(a2)*cos(a3)*sin(a1)];
YXY(2, :) = [sin(a2)*sin(a3), cos(a2), -cos(a3)*sin(a2)];
YXY(3, :) = [- cos(a3)*sin(a1) - cos(a1)*cos(a2)*sin(a3), cos(a1)*sin(a2), cos(a1)*cos(a2)*cos(a3) - sin(a1)*sin(a3)];
dYXYda = simplify(jacobian(YXY(:), [a1 a2 a3]));

%% --------------------------------------------------------------------------------------------

YZY(1, :) = [cos(a1)*cos(a2)*cos(a3) - sin(a1)*sin(a3), -cos(a1)*sin(a2), cos(a3)*sin(a1) + cos(a1)*cos(a2)*sin(a3)];
YZY(2, :) = [cos(a3)*sin(a2), cos(a2), sin(a2)*sin(a3)];
YZY(3, :) = [- cos(a1)*sin(a3) - cos(a2)*cos(a3)*sin(a1), sin(a1)*sin(a2), cos(a1)*cos(a3) - cos(a2)*sin(a1)*sin(a3)];
dYZYda = simplify(jacobian(YZY(:), [a1 a2 a3]));


%% --------------------------------------------------------------------------------------------

ZXZ(1, :) = [cos(a1)*cos(a3) - cos(a2)*sin(a1)*sin(a3), - cos(a1)*sin(a3) - cos(a2)*cos(a3)*sin(a1), sin(a1)*sin(a2)];
ZXZ(2, :) = [cos(a3)*sin(a1) + cos(a1)*cos(a2)*sin(a3), cos(a1)*cos(a2)*cos(a3) - sin(a1)*sin(a3), -cos(a1)*sin(a2)];
ZXZ(3, :) = [sin(a2)*sin(a3), cos(a3)*sin(a2), cos(a2)];
dZXZda = simplify(jacobian(ZXZ(:), [a1 a2 a3]));


%% --------------------------------------------------------------------------------------------

ZYZ(1, :) = [cos(a1)*cos(a2)*cos(a3) - sin(a1)*sin(a3), - cos(a3)*sin(a1) - cos(a1)*cos(a2)*sin(a3), cos(a1)*sin(a2)];
ZYZ(2, :) = [cos(a1)*sin(a3) + cos(a2)*cos(a3)*sin(a1), cos(a1)*cos(a3) - cos(a2)*sin(a1)*sin(a3), sin(a1)*sin(a2)];
ZYZ(3, :) = [-cos(a3)*sin(a2), sin(a2)*sin(a3), cos(a2)];
dZYZda = simplify(jacobian(ZYZ(:), [a1 a2 a3]));


%% --------------------------------------------------------------------------------------------

XYZ(1, :) = [cos(a2)*cos(a3), -cos(a2)*sin(a3), sin(a2)];
XYZ(2, :) = [cos(a1)*sin(a3) + cos(a3)*sin(a1)*sin(a2), cos(a1)*cos(a3) - sin(a1)*sin(a2)*sin(a3), -cos(a2)*sin(a1)];
XYZ(3, :) = [sin(a1)*sin(a3) - cos(a1)*cos(a3)*sin(a2), cos(a3)*sin(a1) + cos(a1)*sin(a2)*sin(a3), cos(a1)*cos(a2)];
dXYZda = simplify(jacobian(XYZ(:), [a1 a2 a3]));


%% --------------------------------------------------------------------------------------------

XZY(1, :) = [cos(a2)*cos(a3), -sin(a2), cos(a2)*sin(a3)];
XZY(2, :) = [sin(a1)*sin(a3) + cos(a1)*cos(a3)*sin(a2), cos(a1)*cos(a2), cos(a1)*sin(a2)*sin(a3) - cos(a3)*sin(a1)];
XZY(3, :) = [cos(a3)*sin(a1)*sin(a2) - cos(a1)*sin(a3), cos(a2)*sin(a1), cos(a1)*cos(a3) + sin(a1)*sin(a2)*sin(a3)];
dXZYda = simplify(jacobian(XZY(:), [a1 a2 a3]));


%% --------------------------------------------------------------------------------------------

YXZ(1, :) = [cos(a1)*cos(a3) + sin(a1)*sin(a2)*sin(a3), cos(a3)*sin(a1)*sin(a2) - cos(a1)*sin(a3), cos(a2)*sin(a1)];
YXZ(2, :) = [cos(a2)*sin(a3), cos(a2)*cos(a3), -sin(a2)];
YXZ(3, :) = [cos(a1)*sin(a2)*sin(a3) - cos(a3)*sin(a1), sin(a1)*sin(a3) + cos(a1)*cos(a3)*sin(a2), cos(a1)*cos(a2)];
dYXZda = simplify(jacobian(YXZ(:), [a1 a2 a3]));


%% --------------------------------------------------------------------------------------------

YZX(1, :) = [cos(a1)*cos(a2), sin(a1)*sin(a3) - cos(a1)*cos(a3)*sin(a2), cos(a3)*sin(a1) + cos(a1)*sin(a2)*sin(a3)];
YZX(2, :) = [sin(a2), cos(a2)*cos(a3), -cos(a2)*sin(a3)];
YZX(3, :) = [-cos(a2)*sin(a1), cos(a1)*sin(a3) + cos(a3)*sin(a1)*sin(a2), cos(a1)*cos(a3) - sin(a1)*sin(a2)*sin(a3)];
dYZXda = simplify(jacobian(YZX(:), [a1 a2 a3]));


%% --------------------------------------------------------------------------------------------

ZXY(1, :) = [cos(a1)*cos(a3) - sin(a1)*sin(a2)*sin(a3), -cos(a2)*sin(a1), cos(a1)*sin(a3) + cos(a3)*sin(a1)*sin(a2)];
ZXY(2, :) = [cos(a3)*sin(a1) + cos(a1)*sin(a2)*sin(a3), cos(a1)*cos(a2), sin(a1)*sin(a3) - cos(a1)*cos(a3)*sin(a2)];
ZXY(3, :) = [-cos(a2)*sin(a3), sin(a2), cos(a2)*cos(a3)];
dZXYda = simplify(jacobian(ZXY(:), [a1 a2 a3]));


%% --------------------------------------------------------------------------------------------

ZYX(1, :) = [cos(a1)*cos(a2), cos(a1)*sin(a2)*sin(a3) - cos(a3)*sin(a1), sin(a1)*sin(a3) + cos(a1)*cos(a3)*sin(a2)];
ZYX(2, :) = [cos(a2)*sin(a1), cos(a1)*cos(a3) + sin(a1)*sin(a2)*sin(a3), cos(a3)*sin(a1)*sin(a2) - cos(a1)*sin(a3)];
ZYX(3, :) = [-sin(a2), cos(a2)*sin(a3), cos(a2)*cos(a3)];
dZYXda = simplify(jacobian(ZYX(:), [a1 a2 a3]));

