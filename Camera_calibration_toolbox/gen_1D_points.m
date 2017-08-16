function [Xp,dXpdXo,dXpdthph,dXpdlen] = gen_1D_points(Xo,thph,rodlen)
% GEN_1D_POINTS compute one dimensional points with origin at Xo, direction angle
% of thph=[theta; phi], and 1D coordinates of rodlen=[x1,x2,...,xn].
% [theta; phi] is the same as the direction of Spherical coordinate system. Theta is
% the angle from Z axis to the vector, and phi is the angle from X axis to the vector;
%
% Input: Xo (3*N or 3*1);     thph (3*N or 3*1)
%
%  Xp = Xo+xn*rodlen;
%
%  See also cartesian2spherical, project_points_mirror2, rigid_refmotion.

% By ZPF @ZVR, 2017-7-19

[m,N] = size(Xo);
assert(ismatrix(Xo) && m==3,'Unexpected dimension of the 1st input!');
[m,n] = size(thph);
assert(ismatrix(thph) && m==2,'Unexpected dimension of the 2nd input!');

sw1 = 0;  % switch to check if Xo have only one column
sw2 = 0;  % switch to check if thph have only one column
if N ~= n,
    if N==1,
        Xo = Xo(:,ones(1,n));
        N = n;
        sw1 = 1;
    elseif n==1,
        thph = thph(:,ones(1,N));
        sw2 = 1;
    else
        error('The columns of the 1st two variables do not match!');
    end;
end;

rodlen = rodlen(:)';
np1D = length(rodlen);
nn = np1D*N;

ctheta = cos(thph(1,:));
stheta = sin(thph(1,:));
cphi = cos(thph(2,:));
sphi = sin(thph(2,:));

% direction of the line.
xn = [stheta.*cphi; stheta.*sphi; ctheta];

Xp = reshape(permute(Xo+xn.*reshape(rodlen,[1,1,np1D]), [1,3,2]), 3,nn);

if nargout > 1,
    flag = N<=5;
    N3 = N*3;
    Nr = np1D*N3;
    n3 = 3*np1D;
    if sw1,
        if flag,
            dXpdXo = repmat(eye(3),nn,1);
        else
            dXpdXo = repmat(speye(3),nn,1);
        end;
    else
        if flag,
            dXpdXo = zeros(Nr,N3);
        else
            dXpdXo = sparse([],[],[],Nr,N3,Nr);
        end;
        for i=1:N,
            dXpdXo((i-1)*n3+1:i*n3, (i-1)*3+1:i*3) = repmat(eye(3),np1D,1);
        end;
    end;
    if nargout > 2,
        one6 = ones(6,1);
        if sw2,
            dxndthph = [ctheta(1)*cphi(1),  -stheta(1)*sphi(1);
                        ctheta(1)*sphi(1),  stheta(1)*cphi(1);
                        -stheta(1),         0];
            dXpdthph = repmat(reshape(rodlen(one6,:),2,n3)'.*repmat(dxndthph,np1D,1),N,1);
        else
            N2 = N*2;
            if flag,
                dXpdthph = zeros(Nr,N2);
            else
                dXpdthph = sparse([],[],[],Nr,N2,Nr*2);
            end;
            for i=1:N,
                dxndthph = [ctheta(i)*cphi(i),  -stheta(i)*sphi(i);
                            ctheta(i)*sphi(i),  stheta(i)*cphi(i);
                            -stheta(i),         0];
                dXpdthph((i-1)*n3+1:i*n3, (i-1)*2+1:i*2) = reshape(rodlen(one6,:),2,n3)'.*repmat(dxndthph,np1D,1);
            end;
        end;
        if nargout > 3,
            if flag,
                dXpdlen = zeros(Nr,np1D);
            else
                dXpdlen = sparse([],[],[],Nr,np1D,Nr);
            end;
            for i=1:np1D,
                j = (i-1)*3;
                dXpdlen(j+1:n3:end,i) = xn(1,:);
                dXpdlen(j+2:n3:end,i) = xn(2,:);
                dXpdlen(j+3:n3:end,i) = xn(3,:);
            end;
        end;
    end;
end;

return;



%% Test of the Jacobians:
m = 5;
np = 1000;
xx = randn(1,m);
tp = [pi*rand(1,np); randn(1,np)];
XX0 = 10*randn(3,np);
[XX,dXXdXX0,dXXdtp,dXXdxx] = gen_1D_points(XX0,tp,xx);
dXX0 = randn(3,np)/100;
dtp = randn(2,np)/100;
dxx = randn(1,m)/100;
XX1 = gen_1D_points(XX0+dXX0,tp+dtp,xx+dxx);
XX_pred = XX + reshape(dXXdXX0*dXX0(:) + dXXdtp*dtp(:) + dXXdxx*dxx(:), 3, []);
gain = norm(XX1-XX)/norm(XX1-XX_pred)

% Test on XX0
XX1 = gen_1D_points(XX0+dXX0,tp,xx);
XX_pred = XX + reshape(dXXdXX0*dXX0(:), 3, []);
gain = norm(XX1-XX)/norm(XX1-XX_pred)

% Test on tp
XX1 = gen_1D_points(XX0,tp+dtp,xx);
XX_pred = XX + reshape(dXXdtp*dtp(:), 3, []);
gain = norm(XX1-XX)/norm(XX1-XX_pred)

% Test on xx
XX1 =gen_1D_points(XX0,tp,xx+dxx);
XX_pred = XX + reshape(dXXdxx*dxx(:), 3, []);
gain = norm(XX1-XX)/norm(XX1-XX_pred)
