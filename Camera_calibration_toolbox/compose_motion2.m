function [OM, T, dOMdOM1,dOMdOM2,dTdOM2,dTdT1,dTdT2] = compose_motion2(OM1,T1,OM2,T2,hand2)
% Given X0 as initial postion.
% X1 = rodrigues(OM1)*diag([1,1,hand1])*X0+T1=R1*h1*X0+T1;
% X2 = rodrigues(OM2)*diag([1,1,hand2])*X1+T2=R2*h2*X1+T2;
% X2 = R2*h2*(R1*h1*X0+T1)+T2=R2*h2*R1*h1*X0+R2*h2*T1+T2
%    = R2*h2*R1*h2*h2*h1*X0+R2*h2*T1+T2
%    = R*h*X0+T
% So R=R2*h2*R1*h2; T=R2*h2*T1+T2; h=h2*h1;
%
% Algorithm:
% Given om=[n1;n2;n3], Hzkk = diag([1,1, -1]), if Rt=Hzkk*rodrigues(om)*Hzkk,
% then we have omt = rodrigues(Rt) = [-n1; -n2; n3]

if nargin<5,
    hand2=1;
end;

if hand2~=1,
    OM1(1:2)=-OM1(1:2);
end;

% hand = hand1*hand2;

if nargout < 3,
    R1 = rodrigues(OM1);       % h2*R1*h2
    R2 = rodrigues(OM2);
    OM = rodrigues(R2 * R1);
    if hand2~=1,
        R2(:,3) = -R2(:,3);
    end;
    T = R2 * T1+ T2;
else
    % Rotations:
    [R1,dR1dOM1] = rodrigues(OM1);
    [R2,dR2dOM2] = rodrigues(OM2);
    R = R2 * R1;
    [dRdR2,dRdR1] = dAB(R2,R1);
    [OM,dOMdR] = rodrigues(R);

    dOMdOM1 = dOMdR * dRdR1 * dR1dOM1;
    dOMdOM2 = dOMdR * dRdR2 * dR2dOM2;
    % Translations:
    if hand2~=1,
        R2(:,3) = -R2(:,3);
    end;
    T = R2 * T1 + T2;
    [dTdR2,dTdT1] = dAB(R2,T1);
    if hand2~=1,
        dOMdOM1(:,1:2) = -dOMdOM1(:,1:2);
        dTdR2(:,7:9) = -dTdR2(:,7:9);
    end;
    dTdOM2 = dTdR2 * dR2dOM2;
    dTdT2 = eye(3);
end;

return;



%% Test of the Jacobians:

om1 = randn(3,1);
om2 = randn(3,1);
T1 = 10*randn(3,1);
T2 = 10*randn(3,1);
hand = sign(randn);
[om3,T3,dom3dom1,dom3dom2,dT3dom2,dT3dT1,dT3dT2] = compose_motion2(om1,T1,om2,T2,hand);

dom1 = randn(3,1) / 1000;
dom2 = randn(3,1) / 1000;
dT1  = randn(3,1) / 10000;
dT2  = randn(3,1) / 10000;

[om3r,T3r] = compose_motion2(om1+dom1,T1+dT1,om2+dom2,T2+dT2,hand);

om3p = om3 + dom3dom1*dom1 +  dom3dT1*dT1 + dom3dom2*dom2 +  dom3dT2*dT2;
T3p  =  T3 + dT3dom1*dom1  +  dT3dT1*dT1  + dT3dom2*dom2  +  dT3dT2*dT2;

norm(om3r - om3) / norm(om3r - om3p)
norm(T3r - T3) / norm(T3r - T3p)

%  Test with inverse composition
om1t = -om1;            % invomw1= rodrigues(Rc1')=-omw1;
if hand~=1,
    om1t(1:2) = -om1t(1:2);      % rodrigues(Hzkk*Rc1'*Hzkk)
end;
Rw1 = rodrigues(om1t);
Rw2 = rodrigues(om2);
omc21 = rodrigues(Rw2*Rw1);
Tc21 = rigid_refmotion(-T1,omc21,T2,hand);

[om21,T21] = compose_motion2(-om1,-rodrigues(-om1)*T1,om2,T2,hand);


