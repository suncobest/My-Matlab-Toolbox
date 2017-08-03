function [Adt, dAdtdom, dAdtdA] = rotspeed2deuler(om, A, op, sw)
% ROTSPEED2DEULER - transform angular velocity to change rate of Euler angle
% depending on Euler angle convention.
%
%  Assuming that Euler angle is defined as intrinsic rotation;
%  OM: treat angular velocity as under intrinsic system (rotating axes), if SW==1;
%    treat angular velocity as under extrinsic system (reference frame), if SW==0.
%  A: Euler angle (orientation);
%  OP: char string of  Euler angle convention;
%        12 options:
%  {'XYX', 'XZX', 'YXY', 'YZY', 'ZXZ', 'ZYZ', 'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
%  SW: mode switch to determine frame system of angular velocity.
%  SW=1 means intrinsic system (default setting); SW=0 means extrinsic system.
%
%  ADT: change rate of Euler angle;
%  DADTDOM: jacobian matrix of ADT over OM.
%  DAdtDA: jacobian matrix of ADT over A.
%   See also deuler2rotspeed, trans_euler_mat, trans_quat_mat, trans_quat_axis, rodrgues.

% by zpf, form BIT, 2016-3-28

[m, n] = size(om);
assert(isequal(size(A),[m,n]) && m==3, 'The 1st and 2nd arguments must be same size with 3 rows!');
if nargin<4 || isempty(sw),
    sw = 1;
else
    sw = ~~sw;
end;

if nargin<3 || isempty(op),
    op = 'ZXZ';
else
    op = upper(op);
    out = double(op);
    assert(ischar(op) && length(out)==3 && all(out>=88) && all(out<=90),'Unexpected input for the 3rd argument!');
end;
n3 = n*3;

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

if n3^2>400,
    dAdtdom = sparse([],[],[],n3,n3,n*9);
else
    dAdtdom = zeros(n3);
end;

if nargout>2,
    % change rate of Euler angle
    om1 = om(1,:);
    om2 = om(2,:);
    om3 = om(3,:);
    dAdtdA = dAdtdom;
end;

switch op,
    case 'XYX',
        if sw,
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [0,      s3(i)/s2(i),      c3(i)/s2(i);
                    0,         c3(i),        -s3(i);
                    1, -c2(i)*s3(i)/s2(i),  -c2(i)*c3(i)/s2(i)];
            end;
        else
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [1, -c2(i)*s1(i)/s2(i),   c1(i)*c2(i)/s2(i);
                    0,       c1(i),           s1(i);
                    0,       s1(i)/s2(i),       -c1(i)/s2(i)];
            end;
        end;
        Adt = reshape(dAdtdom*om(:),3,[]);

        if nargout>2,
            if sw,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [0,   -(c2(i)*(om3(i)*c3(i) + om2(i)*s3(i)))/s2(i)^2,    (om2(i)*c3(i) - om3(i)*s3(i))/s2(i);
                        0,      0,    - om3(i)*c3(i) - om2(i)*s3(i);
                        0,     (om3(i)*c3(i) + om2(i)*s3(i))/s2(i)^2,    -(c2(i)*(om2(i)*c3(i) - om3(i)*s3(i)))/s2(i)];
                end;
            else
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [-(c2(i)*(om2(i)*c1(i) + om3(i)*s1(i)))/s2(i),   (om2(i)*s1(i) - om3(i)*c1(i))/s2(i)^2,   0;
                        om3(i)*c1(i) - om2(i)*s1(i),   0,  0;
                        (om2(i)*c1(i) + om3(i)*s1(i))/s2(i),   (c2(i)*(om3(i)*c1(i) - om2(i)*s1(i)))/s2(i)^2,   0];
                end;
            end;
        end;

    case 'XZX',
        if sw,
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [0,      -c3(i)/s2(i),       s3(i)/s2(i);
                    0,           s3(i),          c3(i);
                    1,   c2(i)*c3(i)/s2(i),  -c2(i)*s3(i)/s2(i)];
            end;
        else
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [1, -c1(i)*c2(i)/s2(i),  -c2(i)*s1(i)/s2(i);
                    0,           -s1(i),          c1(i);
                    0,        c1(i)/s2(i),        s1(i)/s2(i)];
            end;
        end;
        Adt = reshape(dAdtdom*om(:),3,[]);

        if nargout>2,
            if sw,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [0,   (c2(i)*(om2(i)*c3(i) - om3(i)*s3(i)))/s2(i)^2,    (om3(i)*c3(i) + om2(i)*s3(i))/s2(i);
                        0,       0,       om2(i)*c3(i) - om3(i)*s3(i);
                        0,    (om3(i)*s3(i) - om2(i)*c3(i))/s2(i)^2,   -(c2(i)*(om3(i)*c3(i) + om2(i)*s3(i)))/s2(i)];
                end;
            else
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [(c2(i)*(om2(i)*s1(i) - om3(i)*c1(i)))/s2(i),   (om2(i)*c1(i) + om3(i)*s1(i))/s2(i)^2,   0;
                        - om2(i)*c1(i) - om3(i)*s1(i),   0, 0;
                        (om3(i)*c1(i) - om2(i)*s1(i))/s2(i),   -(c2(i)*(om2(i)*c1(i) + om3(i)*s1(i)))/s2(i)^2,   0];
                end;
            end;
        end;

    case 'YXY',
        if sw,
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [s3(i)/s2(i),   0,    -c3(i)/s2(i);
                    c3(i),    0,       s3(i);
                    -c2(i)*s3(i)/s2(i),  1,  c2(i)*c3(i)/s2(i)];
            end;
        else
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [-c2(i)*s1(i)/s2(i), 1,  -c1(i)*c2(i)/s2(i);
                    c1(i),    0,     -s1(i);
                    s1(i)/s2(i),   0,     c1(i)/s2(i)];
            end;
        end;
        Adt = reshape(dAdtdom*om(:),3,[]);

        if nargout>2,
            if sw,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [0,   (c2(i)*(om3(i)*c3(i) - om1(i)*s3(i)))/s2(i)^2,   (om1(i)*c3(i) + om3(i)*s3(i))/s2(i);
                        0,    0,     om3(i)*c3(i) - om1(i)*s3(i);
                        0,    (om1(i)*s3(i) - om3(i)*c3(i))/s2(i)^2,   -(c2(i)*(om1(i)*c3(i) + om3(i)*s3(i)))/s2(i)];
                end;
            else
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [-(c2(i)*(om1(i)*c1(i) - om3(i)*s1(i)))/s2(i),   (om3(i)*c1(i) + om1(i)*s1(i))/s2(i)^2,  0;
                        - om3(i)*c1(i) - om1(i)*s1(i),   0, 0;
                        (om1(i)*c1(i) - om3(i)*s1(i))/s2(i),   -(c2(i)*(om3(i)*c1(i) + om1(i)*s1(i)))/s2(i)^2,   0];
                end;
            end;
        end;

    case 'YZY',
        if sw,
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [c3(i)/s2(i),    0,    s3(i)/s2(i);
                    -s3(i),    0,    c3(i);
                    -c2(i)*c3(i)/s2(i),  1,  -c2(i)*s3(i)/s2(i)];
            end;
        else
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [c1(i)*c2(i)/s2(i),   1,   -c2(i)*s1(i)/s2(i);
                    s1(i),       0,      c1(i);
                    -c1(i)/s2(i),     0,      s1(i)/s2(i)];
            end;
        end;
        Adt = reshape(dAdtdom*om(:),3,[]);

        if nargout>2,
            if sw,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [0,  -(c2(i)*(om1(i)*c3(i) + om3(i)*s3(i)))/s2(i)^2,    (om3(i)*c3(i) - om1(i)*s3(i))/s2(i);
                        0,   0,    - om1(i)*c3(i) - om3(i)*s3(i);
                        0,   (om1(i)*c3(i) + om3(i)*s3(i))/s2(i)^2,   -(c2(i)*(om3(i)*c3(i) - om1(i)*s3(i)))/s2(i)];
                end;
            else
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [-(c2(i)*(om3(i)*c1(i) + om1(i)*s1(i)))/s2(i),   (om3(i)*s1(i) - om1(i)*c1(i))/s2(i)^2,   0;
                        om1(i)*c1(i) - om3(i)*s1(i),   0, 0;
                        (om3(i)*c1(i) + om1(i)*s1(i))/s2(i),    (c2(i)*(om1(i)*c1(i) - om3(i)*s1(i)))/s2(i)^2,  0];
                end;
            end;
        end;

    case 'ZXZ',
        if sw,
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [s3(i)/s2(i),    c3(i)/s2(i),    0;
                    c3(i),     -s3(i),     0;
                    -c2(i)*s3(i)/s2(i), -c2(i)*c3(i)/s2(i), 1];
            end;
        else
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [-c2(i)*s1(i)/s2(i),  c1(i)*c2(i)/s2(i), 1;
                    c1(i),     s1(i),    0;
                    s1(i)/s2(i),     -c1(i)/s2(i),    0];
            end;
        end;
        Adt = reshape(dAdtdom*om(:),3,[]);

        if nargout>2,
            if sw,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [0,   -(c2(i)*(om2(i)*c3(i) + om1(i)*s3(i)))/s2(i)^2,   (om1(i)*c3(i) - om2(i)*s3(i))/s2(i);
                        0,   0,    - om2(i)*c3(i) - om1(i)*s3(i);
                        0,   (om2(i)*c3(i) + om1(i)*s3(i))/s2(i)^2,   -(c2(i)*(om1(i)*c3(i) - om2(i)*s3(i)))/s2(i)];
                end;
            else
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [-(c2(i)*(om1(i)*c1(i) + om2(i)*s1(i)))/s2(i),    (om1(i)*s1(i) - om2(i)*c1(i))/s2(i)^2,  0;
                        om2(i)*c1(i) - om1(i)*s1(i),   0, 0;
                        (om1(i)*c1(i) + om2(i)*s1(i))/s2(i),   (c2(i)*(om2(i)*c1(i) - om1(i)*s1(i)))/s2(i)^2,  0];
                end;
            end;
        end;

    case 'ZYZ',
        if sw,
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [-c3(i)/s2(i),     s3(i)/s2(i),    0;
                    s3(i),     c3(i),     0;
                    c2(i)*c3(i)/s2(i),  -c2(i)*s3(i)/s2(i), 1];
            end;
        else
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [-c1(i)*c2(i)/s2(i),   -c2(i)*s1(i)/s2(i),   1;
                    -s1(i),    c1(i),    0;
                    c1(i)/s2(i),      s1(i)/s2(i),   0];
            end;
        end;
        Adt = reshape(dAdtdom*om(:),3,[]);

        if nargout>2,
            if sw,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [0,   (c2(i)*(om1(i)*c3(i) - om2(i)*s3(i)))/s2(i)^2,   (om2(i)*c3(i) + om1(i)*s3(i))/s2(i);
                        0,    0,     om1(i)*c3(i) - om2(i)*s3(i);
                        0,    (om2(i)*s3(i) - om1(i)*c3(i))/s2(i)^2,  -(c2(i)*(om2(i)*c3(i) + om1(i)*s3(i)))/s2(i)];
                end;
            else
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [-(c2(i)*(om2(i)*c1(i) - om1(i)*s1(i)))/s2(i),   (om1(i)*c1(i) + om2(i)*s1(i))/s2(i)^2,   0;
                        - om1(i)*c1(i) - om2(i)*s1(i),    0,   0;
                        (om2(i)*c1(i) - om1(i)*s1(i))/s2(i),   -(c2(i)*(om1(i)*c1(i) + om2(i)*s1(i)))/s2(i)^2,   0];
                end;
            end;
        end;

    case 'XYZ',
        if sw,
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [c3(i)/c2(i),    -s3(i)/c2(i),    0;
                    s3(i),      c3(i),     0;
                    -c3(i)*s2(i)/c2(i),  s2(i)*s3(i)/c2(i), 1];
            end;
        else
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [1,  s1(i)*s2(i)/c2(i),  -c1(i)*s2(i)/c2(i);
                    0,     c1(i),           s1(i);
                    0,     -s1(i)/c2(i),      c1(i)/c2(i)];
            end;
        end;
        Adt = reshape(dAdtdom*om(:),3,[]);

        if nargout>2,
            if sw,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [0,   (s2(i)*(om1(i)*c3(i) - om2(i)*s3(i)))/c2(i)^2,  -(om2(i)*c3(i) + om1(i)*s3(i))/c2(i);
                        0,    0,    om1(i)*c3(i) - om2(i)*s3(i);
                        0,   -(om1(i)*c3(i) - om2(i)*s3(i))/c2(i)^2,   (s2(i)*(om2(i)*c3(i) + om1(i)*s3(i)))/c2(i)];
                end;
            else
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [(s2(i)*(om2(i)*c1(i) + om3(i)*s1(i)))/c2(i),   -(om3(i)*c1(i) - om2(i)*s1(i))/c2(i)^2,   0;
                        om3(i)*c1(i) - om2(i)*s1(i),   0,  0;
                        -(om2(i)*c1(i) + om3(i)*s1(i))/c2(i),   (s2(i)*(om3(i)*c1(i) - om2(i)*s1(i)))/c2(i)^2,  0];
                end;
            end;
        end;

    case 'XZY',
        if sw,
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [c3(i)/c2(i),   0,   s3(i)/c2(i);
                    -s3(i),     0,     c3(i);
                    c3(i)*s2(i)/c2(i), 1,  s2(i)*s3(i)/c2(i)];
            end;
        else
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [1,   c1(i)*s2(i)/c2(i),   s1(i)*s2(i)/c2(i);
                    0,      -s1(i),         c1(i);
                    0,       c1(i)/c2(i),        s1(i)/c2(i)];
            end;
        end;
        Adt = reshape(dAdtdom*om(:),3,[]);

        if nargout>2,
            if sw,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [0,   (s2(i)*(om1(i)*c3(i) + om3(i)*s3(i)))/c2(i)^2,   (om3(i)*c3(i) - om1(i)*s3(i))/c2(i);
                        0,     0,    - om1(i)*c3(i) - om3(i)*s3(i);
                        0,     (om1(i)*c3(i) + om3(i)*s3(i))/c2(i)^2,   (s2(i)*(om3(i)*c3(i) - om1(i)*s3(i)))/c2(i)];
                end;
            else
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [(s2(i)*(om3(i)*c1(i) - om2(i)*s1(i)))/c2(i),    (om2(i)*c1(i) + om3(i)*s1(i))/c2(i)^2,  0;
                        - om2(i)*c1(i) - om3(i)*s1(i),    0,    0;
                        (om3(i)*c1(i) - om2(i)*s1(i))/c2(i),   (s2(i)*(om2(i)*c1(i) + om3(i)*s1(i)))/c2(i)^2,  0];
                end;
            end;
        end;

    case 'YXZ',
        if sw,
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [s3(i)/c2(i),    c3(i)/c2(i),   0;
                    c3(i),     -s3(i),    0;
                    s2(i)*s3(i)/c2(i),   c3(i)*s2(i)/c2(i),  1];
            end;
        else
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [s1(i)*s2(i)/c2(i),  1,  c1(i)*s2(i)/c2(i);
                    c1(i),     0,     -s1(i);
                    s1(i)/c2(i),     0,     c1(i)/c2(i)];
            end;
        end;
        Adt = reshape(dAdtdom*om(:),3,[]);

        if nargout>2,
            if sw,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [0,  (s2(i)*(om2(i)*c3(i) + om1(i)*s3(i)))/c2(i)^2,    (om1(i)*c3(i) - om2(i)*s3(i))/c2(i);
                        0,     0,     - om2(i)*c3(i) - om1(i)*s3(i);
                        0,     (om2(i)*c3(i) + om1(i)*s3(i))/c2(i)^2,   (s2(i)*(om1(i)*c3(i) - om2(i)*s3(i)))/c2(i)];
                end;
            else
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [(s2(i)*(om1(i)*c1(i) - om3(i)*s1(i)))/c2(i),   (om3(i)*c1(i) + om1(i)*s1(i))/c2(i)^2,  0;
                        - om3(i)*c1(i) - om1(i)*s1(i),   0,  0;
                        (om1(i)*c1(i) - om3(i)*s1(i))/c2(i),   (s2(i)*(om3(i)*c1(i) + om1(i)*s1(i)))/c2(i)^2,  0];
                end;
            end;
        end;

    case 'YZX',
        if sw,
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [0,      c3(i)/c2(i),    -s3(i)/c2(i);
                    0,          s3(i),          c3(i);
                    1,  -c3(i)*s2(i)/c2(i),  s2(i)*s3(i)/c2(i)];
            end;
        else
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [-c1(i)*s2(i)/c2(i), 1,  s1(i)*s2(i)/c2(i);
                    s1(i),      0,       c1(i);
                    c1(i)/c2(i),       0,       -s1(i)/c2(i)];
            end;
        end;
        Adt = reshape(dAdtdom*om(:),3,[]);

        if nargout>2,
            if sw,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [0,   (s2(i)*(om2(i)*c3(i) - om3(i)*s3(i)))/c2(i)^2,   -(om3(i)*c3(i) + om2(i)*s3(i))/c2(i);
                        0,    0,      om2(i)*c3(i) - om3(i)*s3(i);
                        0,   -(om2(i)*c3(i) - om3(i)*s3(i))/c2(i)^2,   (s2(i)*(om3(i)*c3(i) + om2(i)*s3(i)))/c2(i)];
                end;
            else
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [(s2(i)*(om3(i)*c1(i) + om1(i)*s1(i)))/c2(i),   -(om1(i)*c1(i) - om3(i)*s1(i))/c2(i)^2,   0;
                        om1(i)*c1(i) - om3(i)*s1(i),    0,   0;
                        -(om3(i)*c1(i) + om1(i)*s1(i))/c2(i),   (s2(i)*(om1(i)*c1(i) - om3(i)*s1(i)))/c2(i)^2,  0];
                end;
            end;
        end;

    case 'ZXY',
        if sw,
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [-s3(i)/c2(i),     0,      c3(i)/c2(i);
                    c3(i),      0,        s3(i);
                    s2(i)*s3(i)/c2(i), 1, -c3(i)*s2(i)/c2(i)];
            end;
        else
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [s1(i)*s2(i)/c2(i), -c1(i)*s2(i)/c2(i), 1;
                    c1(i),       s1(i),     0;
                    -s1(i)/c2(i),     c1(i)/c2(i),    0];
            end;
        end;
        Adt = reshape(dAdtdom*om(:),3,[]);

        if nargout>2,
            if sw,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [0,   (s2(i)*(om3(i)*c3(i) - om1(i)*s3(i)))/c2(i)^2,   -(om1(i)*c3(i) + om3(i)*s3(i))/c2(i);
                        0,    0,     om3(i)*c3(i) - om1(i)*s3(i);
                        0,    -(om3(i)*c3(i) - om1(i)*s3(i))/c2(i)^2,   (s2(i)*(om1(i)*c3(i) + om3(i)*s3(i)))/c2(i)];
                end;
            else
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [(s2(i)*(om1(i)*c1(i) + om2(i)*s1(i)))/c2(i),   -(om2(i)*c1(i) - om1(i)*s1(i))/c2(i)^2,   0;
                        om2(i)*c1(i) - om1(i)*s1(i),    0,     0;
                        -(om1(i)*c1(i) + om2(i)*s1(i))/c2(i),   (s2(i)*(om2(i)*c1(i) - om1(i)*s1(i)))/c2(i)^2,  0];
                end;
            end;
        end;

    case 'ZYX',
        if sw,
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [0,      s3(i)/c2(i),     c3(i)/c2(i);
                    0,        c3(i),         -s3(i);
                    1,  s2(i)*s3(i)/c2(i),  c3(i)*s2(i)/c2(i)];
            end;
        else
            for i=1:n,
                temp = (i-1)*3;
                ii = temp+1;
                jj = temp+3;
                dAdtdom(ii:jj, ii:jj) = [c1(i)*s2(i)/c2(i),  s1(i)*s2(i)/c2(i), 1;
                    -s1(i),       c1(i),      0;
                    c1(i)/c2(i),     s1(i)/c2(i),     0];
            end;
        end;
        Adt = reshape(dAdtdom*om(:),3,[]);

        if nargout>2,
            if sw,
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [0,  (s2(i)*(om3(i)*c3(i) + om2(i)*s3(i)))/c2(i)^2,   (om2(i)*c3(i) - om3(i)*s3(i))/c2(i);
                        0,     0,      - om3(i)*c3(i) - om2(i)*s3(i);
                        0,     (om3(i)*c3(i) + om2(i)*s3(i))/c2(i)^2, (s2(i)*(om2(i)*c3(i) - om3(i)*s3(i)))/c2(i)];
                end;
            else
                for i=1:n,
                    temp = (i-1)*3;
                    ii = temp+1;
                    jj = temp+3;
                    dAdtdA(ii:jj, ii:jj) = [(s2(i)*(om2(i)*c1(i) - om1(i)*s1(i)))/c2(i),   (om1(i)*c1(i) + om2(i)*s1(i))/c2(i)^2,  0;
                        - om1(i)*c1(i) - om2(i)*s1(i),    0,   0;
                        (om2(i)*c1(i) - om1(i)*s1(i))/c2(i),  (s2(i)*(om1(i)*c1(i) + om2(i)*s1(i)))/c2(i)^2,  0];
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
Rlist = zeros(9,nn);
Rdt = Rlist;
for i=1:nn,
    Rlist(:,i) = reshape(trans_quat_mat(Q1(:,i)),9,[]);
end;
Rdt(:,1) = (Rlist(:,2)-Rlist(:,1))/dt;
Rdt(:,end) = (Rlist(:,end)-Rlist(:,end-1))/dt;
Rdt(:,2:end-1) = (Rlist(:,3:end)-Rlist(:,1:end-2))/(2*dt);
st = {'XYX', 'XZX', 'YXY', 'YZY', 'ZXZ', 'ZYZ', 'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'};
% angular velocity
av0 = zeros(3,nn);
av1 = av0;
for i=1:nn,
    R = reshape(Rlist(:,i),3,[]);
    temp = reshape(Rdt(:,i),3,[])*R';
    temp = (temp-temp')/2;
    av0(:,i) = [temp(3,2); temp(1,3); temp(2,1)];
    av1(:,i) = R'*av0(:,i);
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
    % recompute change rate of Euler angle
    [at0, dat0dav, dat0daa] = rotspeed2deuler(av0, aa, st{k},0);
    err0 = at0-at;
    [at1, dat1dav, dat1daa] = rotspeed2deuler(av1, aa, st{k});
    err1 = at1-at;
    fprintf(1,['\nResult of ' st{k} ' :\n\nnorm(at0-at)=%f;   norm(at1-at)=%f;\n'],[norm(err0),norm(err1)]);

    % test jacobian
    dav = av0/1000;
    daa = aa/1000;
    at01 = rotspeed2deuler(av0+dav, aa+daa, st{k},0);
    at02 = at0+reshape(dat0dav*dav(:),3,[])+reshape(dat0daa*daa(:),3,[]);
    gain0 = norm(at01-at0)/norm(at01-at02);

    at11 = rotspeed2deuler(av1+dav, aa+daa, st{k});
    at12 = at1+reshape(dat1dav*dav(:),3,[])+reshape(dat1daa*daa(:),3,[]);
    gain1 = norm(at11-at1)/norm(at11-at12);
    fprintf(1,'\ngain1=%f;   gain0=%f;\n',[gain1,gain0]);
end;
