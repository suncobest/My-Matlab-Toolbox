function [om3,T3,hand3,dom3dom1,dom3dT1,dom3dom2,dom3dT2,dT3dom1,dT3dT1,dT3dom2,dT3dT2] = compose_transformation(om1,T1,om2,T2,hand1,hand2)

% rotation and inversion

switch nargin
    case 4
        hand1 = 1;
        hand2 = 1;
    case 5
       error('Not eough arguments!')
end

% Rotations:

[R1,dR1dom1] = rodrigues(om1);
R1 = hand1 * R1;
dR1dom1 = hand1 * dR1dom1;

[R2,dR2dom2] =  rodrigues(om2);
R2 = hand2 * R2;
dR2dom2 = hand2 * dR2dom2;

R3 = R2 * R1;

[dR3dR2,dR3dR1] = dAB(R2,R1);

hand3 =  hand1 * hand2;
[om3,dom3dR3] = rodrigues(hand3 * R3);
dom3dR3 = hand3 * dom3dR3;

dom3dom1 = dom3dR3 * dR3dR1 * dR1dom1;
dom3dom2 = dom3dR3 * dR3dR2 * dR2dom2;

dom3dT1 = zeros(3,3);
dom3dT2 = zeros(3,3);


% Translations:

T3t = R2 * T1;

[dT3tdR2,dT3tdT1] = dAB(R2,T1);

dT3tdom2 = dT3tdR2 * dR2dom2;


T3 = T3t + T2;

dT3dT1 = dT3tdT1;
dT3dT2 = eye(3);

dT3dom2 = dT3tdom2;
dT3dom1 = zeros(3,3);


return;

% Test of the Jacobians:

om1 = randn(3,1);
om2 = randn(3,1);
T1 = 10*randn(3,1);
T2 = 10*randn(3,1);
hand1 = sign(randn(1));
hand2 = sign(randn(1));

[om3,T3,hand3,dom3dom1,dom3dT1,dom3dom2,dom3dT2,dT3dom1,dT3dT1,dT3dom2,dT3dT2] = compose_transformation(om1,T1,om2,T2,hand1,hand2);

dom1 = randn(3,1) / 1000;
dom2 = randn(3,1) / 1000;
dT1  = randn(3,1) / 10000;
dT2  = randn(3,1) / 10000;

[om3r,T3r] = compose_transformation(om1+dom1,T1+dT1,om2+dom2,T2+dT2,hand1,hand2);

om3p = om3 + dom3dom1*dom1 +  dom3dT1*dT1 + dom3dom2*dom2 +  dom3dT2*dT2;
T3p  =  T3 + dT3dom1*dom1  +  dT3dT1*dT1  + dT3dom2*dom2  +  dT3dT2*dT2;

norm(om3r - om3) / norm(om3r - om3p)   % 分子是函数变化量，分母是函数余项（二阶以上的高阶项）。大量/小量
norm(T3r - T3) / norm(T3r - T3p)
