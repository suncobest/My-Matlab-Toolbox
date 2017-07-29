function [om1, om0, dom1dAdt, dom1dA, dom0dAdt, dom0dA] = deuler2rotspeed(Adt, A, op)
% DEULER2ROTSPEED - transform change rate of Euler angle to angular velocity
% depending on Euler angle convention.
%
%  Assuming that Euler angle is defined as intrinsic rotation;
%  ADT: change rate of Euler angle;
%  A: Euler angle (orientation);
%  OP: char string of  Euler angle convention;
%        12 options:
%  {'XYX', 'XZX', 'YXY', 'YZY', 'ZXZ', 'ZYZ', 'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
%
%  OM1: angular velocity under intrinsic system (rotating axes)
%  OM0: angular velocity under extrinsic system (reference frame)
%  DOM1DADT: jacobian matrix of OM1 over ADT.
%  DOM1DA: jacobian matrix of OM1 over A.
%  DOM0DADT: jacobian matrix of OM0 over ADT.
%  DOM0DA: jacobian matrix of OM0 over A.
%   See also rotspeed2deuler, trans_euler_mat, trans_quat_mat, trans_quat_axis, rodrgues.

% by zpf, form BIT, 2016-3-27

[m, n] = size(Adt);
assert(isequal(size(A),[m,n]) && m==3, 'The 1st and 2nd arguments must be same size with 3 rows!');
if nargin<3 || isempty(op),
    op = 'ZXZ';
else
    op = upper(op);
    out = double(op);
    assert(ischar(op) && length(out)==3 && all(out>=88) && all(out<=90),'Unexpected input for the 3rd argument!');
end;
n3 = n*3;

% change rate of Euler angle
a1t = Adt(1,:);
a2t = Adt(2,:);
a3t = Adt(3,:);
% Euler angle
a1 = A(1,:);
a2 = A(2,:);
a3 = A(3,:);
% compute sin and cos
s1 = sin(a1);
s2 = sin(a2);
s3 = sin(a3);
c1 = cos(a1);
c2 = cos(a2);
c3 = cos(a3);

if nargout >2,
    if n3^2>400,
        dom1dAdt = sparse([],[],[],n3,n3,n*9);
    else
        dom1dAdt = zeros(n3);
    end;
    dom1dA = dom1dAdt;
    dom0dAdt = dom1dAdt;
    dom0dA = dom1dAdt;
end;

switch op,
    case 'XYX',
        om1 = [a3t + a1t.*c2;
            a2t.*c3 + a1t.*s2.*s3;
            a1t.*c3.*s2 - a2t.*s3];
        
        if nargout>1,
            om0 =  [a1t + a3t.*c2;
                a2t.*c1 + a3t.*s1.*s2;
                a2t.*s1 - a3t.*c1.*s2];
            
            if nargout >2,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    
                    dom1dAdt(ii:jj, ii:jj) = [c2(i),   0, 1;
                        s2(i)*s3(i),  c3(i), 0;
                        c3(i)*s2(i), -s3(i), 0];
                    
                    dom1dA(ii:jj, ii:jj) = [0,     -a1t(i)*s2(i),        0;
                        0, a1t(i)*c2(i)*s3(i),   a1t(i)*c3(i)*s2(i) - a2t(i)*s3(i);
                        0, a1t(i)*c2(i)*c3(i), - a2t(i)*c3(i) - a1t(i)*s2(i)*s3(i)];
                    
                    dom0dAdt(ii:jj, ii:jj) = [1,     0,     c2(i);
                        0, c1(i),  s1(i)*s2(i);
                        0, s1(i), -c1(i)*s2(i)];
                    
                    dom0dA(ii:jj, ii:jj) = [0,   -a3t(i)*s2(i), 0;
                        a3t(i)*c1(i)*s2(i) - a2t(i)*s1(i),  a3t(i)*c2(i)*s1(i), 0;
                        a2t(i)*c1(i) + a3t(i)*s1(i)*s2(i), -a3t(i)*c1(i)*c2(i), 0];
                end;
            end;
        end;
        
    case 'XZX',
        om1 = [a3t + a1t.*c2;
            a2t.*s3 - a1t.*c3.*s2;
            a2t.*c3 + a1t.*s2.*s3];
        
        if nargout>1,
            om0 = [a1t + a3t.*c2;
                a3t.*c1.*s2 - a2t.*s1;
                a2t.*c1 + a3t.*s1.*s2];
            
            if nargout >2,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dom1dAdt(ii:jj, ii:jj) = [c2(i),    0, 1;
                        -c3(i)*s2(i), s3(i), 0;
                        s2(i)*s3(i), c3(i), 0];
                    
                    dom1dA(ii:jj, ii:jj) = [ 0,   -a1t(i)*s2(i),      0;
                        0, -a1t(i)*c2(i)*c3(i), a2t(i)*c3(i) + a1t(i)*s2(i)*s3(i);
                        0,  a1t(i)*c2(i)*s3(i), a1t(i)*c3(i)*s2(i) - a2t(i)*s3(i)];
                    
                    dom0dAdt(ii:jj, ii:jj) = [ 1,     0,     c2(i);
                        0, -s1(i), c1(i)*s2(i);
                        0,  c1(i), s1(i)*s2(i)];
                    
                    dom0dA(ii:jj, ii:jj) = [0,  -a3t(i)*s2(i), 0;
                        - a2t(i)*c1(i) - a3t(i)*s1(i)*s2(i), a3t(i)*c1(i)*c2(i), 0;
                        a3t(i)*c1(i)*s2(i) - a2t(i)*s1(i), a3t(i)*c2(i)*s1(i), 0];
                end;
            end;
        end;
        
    case 'YXY',
        om1 = [a2t.*c3 + a1t.*s2.*s3;
            a3t + a1t.*c2;
            a2t.*s3 - a1t.*c3.*s2];
        
        if nargout>1,
            om0 = [a2t.*c1 + a3t.*s1.*s2;
                a1t + a3t.*c2;
                a3t.*c1.*s2 - a2t.*s1];
            
            if nargout >2,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dom1dAdt(ii:jj, ii:jj) = [s2(i)*s3(i), c3(i), 0;
                        c2(i),    0, 1;
                        -c3(i)*s2(i), s3(i), 0];
                    
                    dom1dA(ii:jj, ii:jj) = [0,  a1t(i)*c2(i)*s3(i), a1t(i)*c3(i)*s2(i) - a2t(i)*s3(i);
                        0,   -a1t(i)*s2(i),     0;
                        0, -a1t(i)*c2(i)*c3(i), a2t(i)*c3(i) + a1t(i)*s2(i)*s3(i)];
                    
                    dom0dAdt(ii:jj, ii:jj) = [0,  c1(i), s1(i)*s2(i);
                        1,        0,         c2(i);
                        0, -s1(i), c1(i)*s2(i)];
                    
                    dom0dA(ii:jj, ii:jj) = [a3t(i)*c1(i)*s2(i) - a2t(i)*s1(i), a3t(i)*c2(i)*s1(i), 0;
                        0,    -a3t(i)*s2(i), 0;
                        - a2t(i)*c1(i) - a3t(i)*s1(i)*s2(i), a3t(i)*c1(i)*c2(i), 0];
                end;
            end;
        end;
        
    case 'YZY',
        om1 = [a1t.*c3.*s2 - a2t.*s3;
            a3t + a1t.*c2;
            a2t.*c3 + a1t.*s2.*s3];
        
        if nargout>1,
            om0 = [a2t.*s1 - a3t.*c1.*s2;
                a1t + a3t.*c2;
                a2t.*c1 + a3t.*s1.*s2];
            
            if nargout >2,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dom1dAdt(ii:jj, ii:jj) = [c3(i)*s2(i), -s3(i), 0;
                        c2(i),    0, 1;
                        s2(i)*s3(i),  c3(i), 0];
                    
                    dom1dA(ii:jj, ii:jj) = [0, a1t(i)*c2(i)*c3(i), - a2t(i)*c3(i) - a1t(i)*s2(i)*s3(i);
                        0,        -a1t(i)*s2(i),            0;
                        0, a1t(i)*c2(i)*s3(i),   a1t(i)*c3(i)*s2(i) - a2t(i)*s3(i)];
                    
                    dom0dAdt(ii:jj, ii:jj) = [0, s1(i), -c1(i)*s2(i);
                        1,       0,          c2(i);
                        0, c1(i),  s1(i)*s2(i)];
                    
                    dom0dA(ii:jj, ii:jj) = [a2t(i)*c1(i) + a3t(i)*s1(i)*s2(i), -a3t(i)*c1(i)*c2(i), 0;
                        0,    -a3t(i)*s2(i), 0;
                        a3t(i)*c1(i)*s2(i) - a2t(i)*s1(i),  a3t(i)*c2(i)*s1(i), 0];
                end;
            end;
        end;
        
    case 'ZXZ',
        om1 = [a2t.*c3 + a1t.*s2.*s3;
            a1t.*c3.*s2 - a2t.*s3;
            a3t + a1t.*c2];
        
        if nargout>1,
            om0 = [a2t.*c1 + a3t.*s1.*s2;
                a2t.*s1 - a3t.*c1.*s2;
                a1t + a3t.*c2];
            
            if nargout >2,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dom1dAdt(ii:jj, ii:jj) = [s2(i)*s3(i),  c3(i), 0;
                        c3(i)*s2(i), -s3(i), 0;
                        c2(i),   0, 1];
                    
                    dom1dA(ii:jj, ii:jj) = [0, a1t(i)*c2(i)*s3(i),   a1t(i)*c3(i)*s2(i) - a2t(i)*s3(i);
                        0, a1t(i)*c2(i)*c3(i), - a2t(i)*c3(i) - a1t(i)*s2(i)*s3(i);
                        0,        -a1t(i)*s2(i),               0];
                    
                    dom0dAdt(ii:jj, ii:jj) = [0, c1(i),  s1(i)*s2(i);
                        0, s1(i), -c1(i)*s2(i);
                        1,       0,        c2(i)];
                    
                    dom0dA(ii:jj, ii:jj) = [a3t(i)*c1(i)*s2(i) - a2t(i)*s1(i),  a3t(i)*c2(i)*s1(i), 0;
                        a2t(i)*c1(i) + a3t(i)*s1(i)*s2(i), -a3t(i)*c1(i)*c2(i), 0;
                        0,   -a3t(i)*s2(i), 0];
                end;
            end;
        end;
        
    case 'ZYZ',
        om1 = [a2t.*s3 - a1t.*c3.*s2;
            a2t.*c3 + a1t.*s2.*s3;
            a3t + a1t.*c2];
        
        if nargout>1,
            om0 = [a3t.*c1.*s2 - a2t.*s1;
                a2t.*c1 + a3t.*s1.*s2;
                a1t + a3t.*c2];
            
            if nargout >2,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dom1dAdt(ii:jj, ii:jj) = [-c3(i)*s2(i), s3(i), 0;
                        s2(i)*s3(i), c3(i), 0;
                        c2(i),    0, 1];
                    
                    dom1dA(ii:jj, ii:jj) = [0, -a1t(i)*c2(i)*c3(i), a2t(i)*c3(i) + a1t(i)*s2(i)*s3(i);
                        0,  a1t(i)*c2(i)*s3(i), a1t(i)*c3(i)*s2(i) - a2t(i)*s3(i);
                        0,         -a1t(i)*s2(i),              0];
                    
                    dom0dAdt(ii:jj, ii:jj) = [0, -s1(i), c1(i)*s2(i);
                        0,  c1(i), s1(i)*s2(i);
                        1,        0,         c2(i)];
                    
                    dom0dA(ii:jj, ii:jj) = [- a2t(i)*c1(i) - a3t(i)*s1(i)*s2(i), a3t(i)*c1(i)*c2(i), 0;
                        a3t(i)*c1(i)*s2(i) - a2t(i)*s1(i), a3t(i)*c2(i)*s1(i), 0;
                        0,    -a3t(i)*s2(i), 0];
                end;
            end;
        end;
        
    case 'XYZ',
        om1 = [a2t.*s3 + a1t.*c2.*c3;
            a2t.*c3 - a1t.*c2.*s3;
            a3t + a1t.*s2];
        
        if nargout>1,
            om0 = [a1t + a3t.*s2;
                a2t.*c1 - a3t.*c2.*s1;
                a2t.*s1 + a3t.*c1.*c2];
            
            if nargout >2,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dom1dAdt(ii:jj, ii:jj) = [c2(i)*c3(i), s3(i), 0;
                        -c2(i)*s3(i), c3(i), 0;
                        s2(i),    0, 1];
                    
                    dom1dA(ii:jj, ii:jj) = [0, -a1t(i)*c3(i)*s2(i),   a2t(i)*c3(i) - a1t(i)*c2(i)*s3(i);
                        0,  a1t(i)*s2(i)*s3(i), - a2t(i)*s3(i) - a1t(i)*c2(i)*c3(i);
                        0,          a1t(i)*c2(i),            0];
                    
                    dom0dAdt(ii:jj, ii:jj) = [1,       0,       s2(i);
                        0, c1(i), -c2(i)*s1(i);
                        0, s1(i),  c1(i)*c2(i)];
                    
                    dom0dA(ii:jj, ii:jj) = [0,     a3t(i)*c2(i), 0;
                        - a2t(i)*s1(i) - a3t(i)*c1(i)*c2(i),  a3t(i)*s1(i)*s2(i), 0;
                        a2t(i)*c1(i) - a3t(i)*c2(i)*s1(i), -a3t(i)*c1(i)*s2(i), 0];
                end;
            end;
        end;
        
    case 'XZY',
        om1 = [ a1t.*c2.*c3 - a2t.*s3;
            a3t - a1t.*s2;
            a2t.*c3 + a1t.*c2.*s3];
        
        if nargout>1,
            om0 = [a1t - a3t.*s2;
                a3t.*c1.*c2 - a2t.*s1;
                a2t.*c1 + a3t.*c2.*s1];
            
            if nargout >2,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dom1dAdt(ii:jj, ii:jj) = [c2(i)*c3(i), -s3(i), 0;
                        -s2(i),   0, 1;
                        c2(i)*s3(i),  c3(i), 0];
                    
                    dom1dA(ii:jj, ii:jj) = [0, -a1t(i)*c3(i)*s2(i), - a2t(i)*c3(i) - a1t(i)*c2(i)*s3(i);
                        0,         -a1t(i)*c2(i),          0;
                        0, -a1t(i)*s2(i)*s3(i),   a1t(i)*c2(i)*c3(i) - a2t(i)*s3(i)];
                    
                    dom0dAdt(ii:jj, ii:jj) = [1,      0,      -s2(i);
                        0, -s1(i), c1(i)*c2(i);
                        0,  c1(i), c2(i)*s1(i)];
                    
                    dom0dA(ii:jj, ii:jj) = [0,    -a3t(i)*c2(i), 0;
                        - a2t(i)*c1(i) - a3t(i)*c2(i)*s1(i), -a3t(i)*c1(i)*s2(i), 0;
                        a3t(i)*c1(i)*c2(i) - a2t(i)*s1(i), -a3t(i)*s1(i)*s2(i), 0];
                end;
            end;
        end;
        
    case 'YXZ',
        om1 = [a2t.*c3 + a1t.*c2.*s3;
            a1t.*c2.*c3 - a2t.*s3;
            a3t - a1t.*s2];
        
        if nargout>1,
            om0 = [a2t.*c1 + a3t.*c2.*s1;
                a1t - a3t.*s2;
                a3t.*c1.*c2 - a2t.*s1];
            
            if nargout >2,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dom1dAdt(ii:jj, ii:jj) = [c2(i)*s3(i),  c3(i), 0;
                        c2(i)*c3(i), -s3(i), 0;
                        -s2(i),   0, 1];
                    
                    dom1dA(ii:jj, ii:jj) = [0, -a1t(i)*s2(i)*s3(i),   a1t(i)*c2(i)*c3(i) - a2t(i)*s3(i);
                        0, -a1t(i)*c3(i)*s2(i), - a2t(i)*c3(i) - a1t(i)*c2(i)*s3(i);
                        0,       -a1t(i)*c2(i),           0];
                    
                    dom0dAdt(ii:jj, ii:jj) = [0,  c1(i), c2(i)*s1(i);
                        1,        0,       -s2(i);
                        0, -s1(i), c1(i)*c2(i)];
                    
                    dom0dA(ii:jj, ii:jj) = [a3t(i)*c1(i)*c2(i) - a2t(i)*s1(i), -a3t(i)*s1(i)*s2(i), 0;
                        0,   -a3t(i)*c2(i), 0;
                        - a2t(i)*c1(i) - a3t(i)*c2(i)*s1(i), -a3t(i)*c1(i)*s2(i), 0];
                end;
            end;
        end;
        
    case 'YZX',
        om1 = [a3t + a1t.*s2;
            a2t.*s3 + a1t.*c2.*c3;
            a2t.*c3 - a1t.*c2.*s3];
        
        if nargout>1,
            om0 = [a2t.*s1 + a3t.*c1.*c2;
                a1t + a3t.*s2;
                a2t.*c1 - a3t.*c2.*s1];
            
            if nargout >2,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dom1dAdt(ii:jj, ii:jj) = [s2(i),   0, 1;
                        c2(i)*c3(i), s3(i), 0;
                        -c2(i)*s3(i), c3(i), 0];
                    
                    dom1dA(ii:jj, ii:jj) = [0,          a1t(i)*c2(i),             0;
                        0, -a1t(i)*c3(i)*s2(i),   a2t(i)*c3(i) - a1t(i)*c2(i)*s3(i);
                        0,  a1t(i)*s2(i)*s3(i), - a2t(i)*s3(i) - a1t(i)*c2(i)*c3(i)];
                    
                    dom0dAdt(ii:jj, ii:jj) = [0, s1(i),  c1(i)*c2(i);
                        1,       0,         s2(i);
                        0, c1(i), -c2(i)*s1(i)];
                    
                    dom0dA(ii:jj, ii:jj) = [a2t(i)*c1(i) - a3t(i)*c2(i)*s1(i), -a3t(i)*c1(i)*s2(i), 0;
                        0,    a3t(i)*c2(i), 0;
                        - a2t(i)*s1(i) - a3t(i)*c1(i)*c2(i),  a3t(i)*s1(i)*s2(i), 0];
                end;
            end;
        end;
        
    case 'ZXY',
        om1 = [a2t.*c3 - a1t.*c2.*s3;
            a3t + a1t.*s2;
            a2t.*s3 + a1t.*c2.*c3];
        
        if nargout>1,
            om0 = [a2t.*c1 - a3t.*c2.*s1;
                a2t.*s1 + a3t.*c1.*c2;
                a1t + a3t.*s2];
            
            if nargout >2,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dom1dAdt(ii:jj, ii:jj) = [-c2(i)*s3(i), c3(i), 0;
                        s2(i),    0, 1;
                        c2(i)*c3(i), s3(i), 0];
                    
                    dom1dA(ii:jj, ii:jj) = [0,  a1t(i)*s2(i)*s3(i), - a2t(i)*s3(i) - a1t(i)*c2(i)*c3(i);
                        0,          a1t(i)*c2(i),               0;
                        0, -a1t(i)*c3(i)*s2(i),   a2t(i)*c3(i) - a1t(i)*c2(i)*s3(i)];
                    
                    dom0dAdt(ii:jj, ii:jj) = [0, c1(i), -c2(i)*s1(i);
                        0, s1(i),  c1(i)*c2(i);
                        1,       0,        s2(i)];
                    
                    dom0dA(ii:jj, ii:jj) = [- a2t(i)*s1(i) - a3t(i)*c1(i)*c2(i),  a3t(i)*s1(i)*s2(i), 0;
                        a2t(i)*c1(i) - a3t(i)*c2(i)*s1(i), -a3t(i)*c1(i)*s2(i), 0;
                        0,     a3t(i)*c2(i), 0];
                end;
            end;
        end;
        
    case 'ZYX',
        om1 = [a3t - a1t.*s2;
            a2t.*c3 + a1t.*c2.*s3;
            a1t.*c2.*c3 - a2t.*s3];
        
        if nargout>1,
            om0 = [a3t.*c1.*c2 - a2t.*s1;
                a2t.*c1 + a3t.*c2.*s1;
                a1t - a3t.*s2];
            
            if nargout >2,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dom1dAdt(ii:jj, ii:jj) = [-s2(i),   0, 1;
                        c2(i)*s3(i),  c3(i), 0;
                        c2(i)*c3(i), -s3(i), 0];
                    
                    dom1dA(ii:jj, ii:jj) = [0,         -a1t(i)*c2(i),           0;
                        0, -a1t(i)*s2(i)*s3(i),   a1t(i)*c2(i)*c3(i) - a2t(i)*s3(i);
                        0, -a1t(i)*c3(i)*s2(i), - a2t(i)*c3(i) - a1t(i)*c2(i)*s3(i)];
                    
                    dom0dAdt(ii:jj, ii:jj) = [0, -s1(i), c1(i)*c2(i);
                        0,  c1(i), c2(i)*s1(i);
                        1,        0,        -s2(i)];
                    
                    dom0dA(ii:jj, ii:jj) = [- a2t(i)*c1(i) - a3t(i)*c2(i)*s1(i), -a3t(i)*c1(i)*s2(i), 0;
                        a3t(i)*c1(i)*c2(i) - a2t(i)*s1(i), -a3t(i)*s1(i)*s2(i), 0;
                        0,   -a3t(i)*c2(i), 0];
                end;
            end;
        end;
        
    otherwise,
        error('Unexpected input for the 3rd argument!');
end;

return;
        



%% Test of function:
n = 3;
nn = 2001;
om = randn(3,n);
t0 = (0:n-1)/n-1;
t1 = linspace(0,1,nn);
dt = t1(2)-t1(1);
Q0 = trans_quat_axis(om);
Q1 = squad(t0,Q0,t1);
st = {'XYX', 'XZX', 'YXY', 'YZY', 'ZXZ', 'ZYZ', 'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'};
% angular velocity
av = zeros(3,nn);
Rlist = zeros(9,nn);
Rdt = Rlist;
for i=1:nn,
    Rlist(:,i) = reshape(trans_quat_mat(Q1(:,i)),9,[]);
end;
Rdt(:,1) = (Rlist(:,2)-Rlist(:,1))/dt;
Rdt(:,end) = (Rlist(:,end)-Rlist(:,end-1))/dt;
Rdt(:,2:end-1) = (Rlist(:,3:end)-Rlist(:,1:end-2))/(2*dt);
for i=1:nn,
    temp = reshape(Rdt(:,i),3,[])*reshape(Rlist(:,i),3,[])';
    temp = (temp-temp')/2;
    av(:,i) = [temp(3,2); temp(1,3); temp(2,1)];
end;

for k=1:12,
    % Euler angle
    aa = zeros(3,nn);
    at = aa;
    for i=1:nn,
        aa(:,i) = trans_euler_mat(reshape(Rlist(:,i),3,[]),st{k});
    end;
    % change rate of Euler angle
    at(:,1) = aa(:,2)-aa(:,1);
    at(:,end) = aa(:,end)-aa(:,end-1);
    at(:,2:end-1) = aa(:,3:end)-aa(:,1:end-2);
    id = at>pi;
    at(id) = at(id)-2*pi;
    id = at<-pi;
    at(id) = at(id)+2*pi;
    at(:,2:end-1) = at(:,2:end-1)/2;
    at = at/dt;
    % recompute angular velocity
    [av1,av0,dav1dat, dav1daa, dav0dat, dav0daa] = deuler2rotspeed(at, aa, st{k});
    err = av0-av;
    fprintf(1,['\nResult of ' st{k} ' :\n\nnorm(av0-av)=%f\n'],norm(err));
    % test jacobian
    dat = at/1000;
    daa = aa/1000;
    [av11, av01] = deuler2rotspeed(at+dat, aa+daa, st{k});
    av12 = av1+reshape(dav1dat*dat(:),3,[])+reshape(dav1daa*daa(:),3,[]);
    gain1 = norm(av11-av1)/norm(av11-av12);
    av02 = av0+reshape(dav0dat*dat(:),3,[])+reshape(dav0daa*daa(:),3,[]);
    gain0 = norm(av01-av0)/norm(av01-av02);
    fprintf(1,'\ngain1=%f;   gain0=%f;\n',[gain1,gain0]);
end;
