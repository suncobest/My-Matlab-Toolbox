function [Q, dQdA, dQdB] = quatmul(A, B, SW)
% QUATMUL - Computes product of two groups of quaternions
%
% Usage: Q = quatmul(A, B), multiply B by A on the left side.
%
% Arguments: A, B - Quaternions assumed to be 4-vectors in the
%                   form  A = A(1:4,1:na), B = B(1:4,1:nb); SW is the
%                   switch to turn on the function of multiplying
%                   corresponding quaternion element if SW==1.
% Returns:   Q =  Q(1:4, 1:nb, 1:na);   - Quaternion products in the form
%        if nb==1,  Q(:,i) = quatmul(A(:,i),B) ; if nb>1, Q(:,:,i) = quatmul(A(:,i),B)
%
%   dQdA is the derivative of Q wrt A, dQdB is the derivative of Q wrt B
%   wrt : with respect to
%
% If Q is the product of quaternion A and quaternion B, then Q = AB
% Q(1)  =   A(4)*B(1) - A(3)*B(2) + A(2)*B(3) + A(1)*B(4);
% Q(2)  =   A(3)*B(1) + A(4)*B(2) - A(1)*B(3) + A(2)*B(4);
% Q(3)  =-A(2)*B(1) + A(1)*B(2) + A(4)*B(3) + A(3)*B(4);
% Q(4)  =-A(1)*B(1) - A(2)*B(2) - A(3)*B(3) + A(4)*B(4);
% See also: trans_quat_mat, trans_quat_axis, slerp.

assert(ismatrix(A) && ismatrix(B), 'The two arguments (quaternions) must be matrix!');
[ma, na]=size(A);
switch ma,
    case 1,
        switch na,
            case 3,
                A = [A,0]';
                na = 1;
            case 4,
                A = A';
                na = 1;
            otherwise,
                error('Unexpected dimension of the 1st argument!');
        end;
    case 3,
        A = [A; zeros(1, na)];
    case 4;
    otherwise,
        error('Unexpected dimension of the 1st argument!');
end;

[mb, nb]=size(B);
switch mb,
    case 1,
        switch nb,
            case 3,
                B = [B,0]';
                nb = 1;
            case 4,
                B = B';
                nb = 1;
            otherwise,
                error('Unexpected dimension of the 2nd argument!');
        end;
    case 3,
        B = [B; zeros(1, nb)];
    case 4;
    otherwise,
        error('Unexpected dimension of the 2nd argument!');
end;

if nargin<3,
    SW = false;
else
    SW = ~~SW;
end;

if nb ==1,
    B1=B(1);
    B2=B(2);
    B3=B(3);
    B4=B(4);
    % matrix form
    Qb = [B4, B3, -B2, B1;
        -B3, B4, B1, B2;
        B2,-B1, B4, B3;
        -B1, -B2, -B3, B4];
    Q = Qb*A;
    
    if nargout <2,
        return;
    end;   
    ndQ =4*na;    
    if ndQ^2>400,
        dQdA = sparse([],[],[],ndQ,ndQ,16*na);
    else
        dQdA = zeros(ndQ);
    end;
    dQdB = zeros(ndQ,4);
    for i = 1:na,
        ni = 4*(i-1)+1:4*i;
        dQdA(ni,ni) = Qb;
        A1=A(1,i);
        A2=A(2,i);
        A3=A(3,i);
        A4=A(4,i);
        % matrix form
        Qa = [A4, -A3, A2, A1;
            A3, A4, -A1, A2;
            -A2, A1, A4, A3;
            -A1, -A2, -A3, A4];
        dQdB(4*(i-1)+1:4*i, :) = Qa;
    end;

elseif na == 1,
    A1=A(1);
    A2=A(2);
    A3=A(3);
    A4=A(4);
    % matrix form
    Qa = [A4, -A3, A2, A1;
        A3, A4, -A1, A2;
        -A2, A1, A4, A3;
        -A1, -A2, -A3, A4];
    Q = Qa*B;
    
    if nargout <2,
        return;
    end;
    ndQ =4*nb;
    dQdA =  zeros(ndQ,4);
    if ndQ^2>400,
        dQdB = sparse([],[],[],ndQ,ndQ,16*nb);
    else
        dQdB = zeros(ndQ);
    end;
    for i =1:nb,
        ni = 4*(i-1)+1:4*i;
        dQdB(ni,ni) = Qa;
        B1=B(1,i);
        B2=B(2,i);
        B3=B(3,i);
        B4=B(4,i);
        Qb = [B4, B3, -B2, B1;
            -B3, B4, B1, B2;
            B2,-B1, B4, B3;
            -B1, -B2, -B3, B4];
        dQdA(ni, :) = Qb;
    end;
       
elseif na == nb && SW,
    Q = zeros(4,nb);
    ndQ =4*na;
    if ndQ^2>400,
        dQdA = sparse([],[],[],ndQ,ndQ,16*na);
    else
        dQdA = zeros(ndQ);
    end;
    for i = 1:na,
        B1=B(1,i);
        B2=B(2,i);
        B3=B(3,i);
        B4=B(4,i);
        ni = 4*(i-1)+1:4*i;
        % matrix form
        Qb = [B4, B3, -B2, B1;
            -B3, B4, B1, B2;
            B2,-B1, B4, B3;
            -B1, -B2, -B3, B4];
        Q(:,i) = Qb*A(:,i);
        dQdA(ni,ni) = Qb;
    end;
    
    if nargout <3,
        return;
    end;
    if ndQ^2>400,
        dQdB = sparse([],[],[],ndQ,ndQ,16*na);
    else
        dQdB = zeros(ndQ);
    end;
    for i = 1:na,
        A1=A(1,i);
        A2=A(2,i);
        A3=A(3,i);
        A4=A(4,i);
        % matrix form
        Qa = [A4, -A3, A2, A1;
            A3, A4, -A1, A2;
            -A2, A1, A4, A3;
            -A1, -A2, -A3, A4];
        ni = 4*(i-1)+1:4*i;
        dQdB(ni, ni) = Qa;
    end;
    
else
    Q = zeros(4,nb,na);
    for i = 1:na,
        A1=A(1,i);
        A2=A(2,i);
        A3=A(3,i);
        A4=A(4,i);
        % matrix form
        Qi = [A4, -A3, A2, A1;
            A3, A4, -A1, A2;
            -A2, A1, A4, A3;
            -A1, -A2, -A3, A4];
        Q(:,:,i) = Qi*B;
    end;
    if nargout >1,
        fprintf(1,'Jacobian Matrix is not calculated!\n');
        dQdA = [];
        dQdB = [];
        return;
    end;
end;

return;


%% Test
a = randn(4,1);
a = a/norm(a);
inva = [-a(1:3); a(4)];
R = trans_quat_mat(a);
b = randn(3,7);
ab = quatmul(quatmul(a,b),inva);
err = R*b-ab(1:3,:)
%% Test of dQdA, dQdB
m = 4;
n = 4;
sw = 1;
a=randn(4,m);
b=randn(4,n);
da = randn(4,m)/100;
db = randn(4,n)/100;
[c,dcda,dcdb]=quatmul(a,b,sw);
c1=quatmul(a+da,b,sw);
c2=c+reshape(dcda*da(:),4,[]);
gain = norm(c1-c)/norm(c1-c2)
%%
c1=quatmul(a,b+db,sw);
c2=c+reshape(dcdb*db(:),4,[]);
gain = norm(c1-c)/norm(c1-c2)
%%
c1=quatmul(a+da,b+db,sw);
c2=c+reshape(dcda*da(:),4,[])+reshape(dcdb*db(:),4,[]);
gain = norm(c1-c)/norm(c1-c2)

