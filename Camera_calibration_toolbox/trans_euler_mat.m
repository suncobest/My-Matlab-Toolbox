function [out, dout] = trans_euler_mat( in, op)
% TRANS_EULER_AXIS - transform Euler angle to rotation matrix, or vice versa.
%  EULER2MAT transform Euler angle to rotation matrix, while MAT2EULER
%  transform rotation matrix to Euler angle.
%
%  Assuming that Euler angle is defined as intrinsic rotation;
%  IN: Euler angle or rotation matrix;
%  OP: char string of  Euler angle convention;
%        12 options:
%  {'XYX', 'XZX', 'YXY', 'YZY', 'ZXZ', 'ZYZ', 'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'}
%
%  OUT: rotation matrix or Euler angle depending on input.
%  DOUT: jacobian matrix of OUT over IN. 
%   See also rodrgues, trans_quat_mat, trans_quat_axis.

% by zpf, form BIT, 2016-3-25

bigeps = 1e5;     %1e+16*eps;
if nargin<2 || isempty(op),
    op = 'ZXZ';
else
    op = upper(op);
    out = double(op);
    assert(ischar(op) && length(out)==3 && all(out>=88) && all(out<=90),'Unexpected input for the 2nd argument!');
end;

if isvector(in) && length(in)==3,
    if nargout>1,
        [out,dout] = euler2mat(in,op);
    else
        out = euler2mat(in,op);
    end;
elseif isequal(size(in),[3,3]) && norm(in' * in - eye(3)) < bigeps,
    detR = det(in);
    assert( detR>0 && abs(detR-1)<bigeps, 'The input rotation matrix (3x3) must have positive determinant!');
    if nargout>1,
        [out,dout] = mat2euler(in,op);
    else
        out = mat2euler(in,op);
    end;
else
     error('The 1st argument must be either Euler angle vector or rotation matrix!');
end;

return;


function [R, dRdvec] = euler2mat(vec, op)

a1 = vec(1);
a2 = vec(2);
a3 = vec(3);
% compute sin and cos
s1 = sin(a1);
s2 = sin(a2);
s3 = sin(a3);
c1 = cos(a1);
c2 = cos(a2);
c3 = cos(a3);

switch op,
    case 'XYX',
        R = [c2,                                s2*s3,                                     c3*s2;
            s1*s2,                  c1*c3-c2*s1*s3,                     -c1*s3-c2*c3*s1;
            -c1*s2,             c3*s1+c1*c2*s3,                         c1*c2*c3-s1*s3];
        
        if nargout>1,
            dRdvec = [ 0,                                  -s2,                                     0;
                c1*s2,                                c2*s1,                                    0;
                s1*s2,                              -c1*c2,                                   0;
                0,                                       c2*s3,                                  c3*s2;
                -c3*s1-c1*c2*s3,           s1*s2*s3,                -c1*s3-c2*c3*s1;
                c1*c3-c2*s1*s3,           -c1*s2*s3,                  c1*c2*c3-s1*s3;
                0,                                        c2*c3,                            -s2*s3;
                s1*s3-c1*c2*c3,             c3*s1*s2,             c2*s1*s3-c1*c3;
                -c1*s3-c2*c3*s1,        -c1*c3*s2,          -c3*s1-c1*c2*s3];
        end;
        
    case 'XZX',
        R = [c2,                      -c3*s2,                                   s2*s3;
            c1*s2,                   c1*c2*c3-s1*s3,            -c3*s1-c1*c2*s3;
            s1*s2,                 	c1*s3+c2*c3*s1,               c1*c3-c2*s1*s3];
        
        if nargout>1,
            dRdvec = [0,                                 -s2,                                     0;
                -s1*s2,                                  c1*c2,                                   0;
                c1*s2,                                     c2*s1,                                  0;
                0,                                         -c2*c3,                               s2*s3;
                -c1*s3-c2*c3*s1,          -c1*c3*s2,               -c3*s1-c1*c2*s3;
                c1*c2*c3-s1*s3,            -c3*s1*s2,                  c1*c3-c2*s1*s3;
                0,                                           c2*s3,                               c3*s2;
                c2*s1*s3-c1*c3,               c1*s2*s3,                s1*s3-c1*c2*c3;
                -c3*s1-c1*c2*s3,            s1*s2*s3,                -c1*s3-c2*c3*s1];
        end;

  case 'YXY',
        R = [c1*c3-c2*s1*s3,    	      s1*s2,       	    c1*s3+c2*c3*s1;
            s2*s3,                                       c2,                            -c3*s2;
            -c3*s1-c1*c2*s3,	    c1*s2,	    c1*c2*c3-s1*s3];
        
        if nargout>1,
            dRdvec = [-c3*s1-c1*c2*s3,        s1*s2*s3,        -c1*s3-c2*c3*s1;
                0,                                       c2*s3,                          c3*s2;
                c2*s1*s3-c1*c3,         c1*s2*s3,            s1*s3-c1*c2*c3;
                c1*s2,                                c2*s1,                                0;
                0,                                      -s2,                                    0;
                -s1*s2,                            c1*c2,                                 0;
                c1*c2*c3-s1*s3,            -c3*s1*s2,	    c1*c3-c2*s1*s3;
                0,                                     -c2*c3,                          s2*s3;
                -c1*s3-c2*c3*s1,	  -c1*c3*s2,	    -c3*s1-c1*c2*s3];
        end;
        
    case 'YZY',
        R = [c1*c2*c3-s1*s3,            -c1*s2,                c3*s1+c1*c2*s3;
            c3*s2,                                   c2,                                s2*s3;
            -c1*s3-c2*c3*s1,              s1*s2,                   c1*c3-c2*s1*s3];
        
        if nargout>1,
            dRdvec = [-c1*s3-c2*c3*s1,                -c1*c3*s2,                -c3*s1-c1*c2*s3;
                0,                                              c2*c3,                                    -s2*s3;
                s1*s3-c1*c2*c3,                c3*s1*s2,                     c2*s1*s3-c1*c3;
                s1*s2,                                   -c1*c2,                                          0;
                0,                                            -s2,                                              0;
                c1*s2,                                      c2*s1,                                           0;
                c1*c3-c2*s1*s3,                -c1*s2*s3,                    c1*c2*c3-s1*s3;
                0,                                            c2*s3,                                       c3*s2;
                -c3*s1-c1*c2*s3,                s1*s2*s3,                -c1*s3-c2*c3*s1];
        end;
        
    case 'ZXZ',
        R = [c1*c3-c2*s1*s3,                -c1*s3-c2*c3*s1,                s1*s2;
            c3*s1+c1*c2*s3,                c1*c2*c3-s1*s3,                -c1*s2;
            s2*s3,                                 c3*s2,                                        c2];
        
        if nargout>1,
            dRdvec = [-c3*s1-c1*c2*s3,                s1*s2*s3,                -c1*s3-c2*c3*s1;
                c1*c3-c2*s1*s3,                -c1*s2*s3,                c1*c2*c3-s1*s3;
                0,                                         c2*s3,                             c3*s2;
                s1*s3-c1*c2*c3,                c3*s1*s2,                c2*s1*s3 - c1*c3;
                -c1*s3-c2*c3*s1,                -c1*c3*s2,                -c3*s1-c1*c2*s3;
                0,                                         c2*c3,                              -s2*s3;
                c1*s2,                                  c2*s1,                                        0;
                s1*s2,                                 -c1*c2,                                      0;
                0,                                         -s2,                                           0];
        end;
        
    case 'ZYZ',
        R = [c1*c2*c3-s1*s3,	-c3*s1-c1*c2*s3,	c1*s2;
            c1*s3+c2*c3*s1,        c1*c3-c2*s1*s3,         s1*s2;
            -c3*s2,                               s2*s3,                            c2];
        
        if nargout>1,
            dRdvec = [-c1*s3-c2*c3*s1,                -c1*c3*s2,                -c3*s1-c1*c2*s3;
                c1*c2*c3-s1*s3,                -c3*s1*s2,                c1*c3-c2*s1*s3;
                0,                                       -c2*c3,                              s2*s3;
                c2*s1*s3-c1*c3,                c1*s2*s3,                s1*s3-c1*c2*c3;
                -c3*s1-c1*c2*s3,                s1*s2*s3,                -c1*s3-c2*c3*s1;
                0,                                         c2*s3,                               c3*s2;
                -s1*s2,                               c1*c2,                                   0;
                c1*s2,                                  c2*s1,                                   0;
                0,                                          -s2,                                     0];
        end;
        
    case 'XYZ',
        R = [c2*c3,                              -c2*s3,                              s2;
            c1*s3+c3*s1*s2,         c1*c3-s1*s2*s3,            -c2*s1;
            s1*s3-c1*c3*s2,         c3*s1 + c1*s2*s3,           c1*c2];
        
        if nargout>1,
            dRdvec = [0,                        -c3*s2,                          -c2*s3;
                c1*c3*s2-s1*s3,                c2*c3*s1,                c1*c3-s1*s2*s3;
                c1*s3+c3*s1*s2,                -c1*c2*c3,                c3*s1+c1*s2*s3;
                0,                                        s2*s3,                             -c2*c3;
                -c3*s1-c1*s2*s3,                -c2*s1*s3,                -c1*s3-c3*s1*s2;
                c1*c3-s1*s2*s3,                c1*c2*s3,                c1*c3*s2-s1*s3;
                0,                                          c2,                               0;
                -c1*c2,                             s1*s2,                             0;
                -c2*s1,                         -c1*s2,                               0];
        end;
        
    case 'XZY',
        R = [c2*c3,                              -s2,                                   c2*s3;
            s1*s3+c1*c3*s2,                c1*c2,                c1*s2*s3 - c3*s1;
            c3*s1*s2-c1*s3,                c2*s1,                c1*c3+s1*s2*s3];
        
        if nargout>1,
            dRdvec = [0,                        -c3*s2,                           -c2*s3;
                c1*s3-c3*s1*s2,                c1*c2*c3,                c3*s1-c1*s2*s3;
                s1*s3+c1*c3*s2,                c2*c3*s1,                -c1*c3-s1*s2*s3;
                0,                                        -c2,                                      0;
                -c2*s1,                            -c1*s2,                                   0;
                c1*c2,                              -s1*s2,                                    0;
                0,                                      -s2*s3,                             c2*c3;
                -c1*c3-s1*s2*s3,                c1*c2*s3,                s1*s3+c1*c3*s2;
                c1*s2*s3-c3*s1,                c2*s1*s3,                c3*s1*s2-c1*s3];
        end;
        
    case 'YXZ',
        R = [c1*c3+s1*s2*s3,                c3*s1*s2-c1*s3,                c2*s1;
            c2*s3,                                      c2*c3,                               -s2;
            c1*s2*s3-c3*s1,            s1*s3+c1*c3*s2,                c1*c2];
        
        if nargout>1,
            dRdvec = [c1*s2*s3-c3*s1,                c2*s1*s3,                c3*s1*s2-c1*s3;
                0,                                           -s2*s3,                              c2*c3;
                -c1*c3-s1*s2*s3,                c1*c2*s3,                s1*s3+c1*c3*s2;
                s1*s3+c1*c3*s2,                c2*c3*s1,                -c1*c3-s1*s2*s3;
                0,                                           -c3*s2,                             -c2*s3;
                c1*s3-c3*s1*s2,                c1*c2*c3,                c3*s1-c1*s2*s3;
                c1*c2,                                    -s1*s2,                                      0;
                0,                                              -c2,                                         0;
                -c2*s1,                                  -c1*s2,                                     0];
        end;
        
    case 'YZX',
        R = [c1*c2,                s1*s3-c1*c3*s2,                c3*s1+c1*s2*s3;
            s2,                                  c2*c3,                                     -c2*s3;
            -c2*s1,                c1*s3+c3*s1*s2,                c1*c3-s1*s2*s3];
        
        if nargout>1,
            dRdvec = [-c2*s1,                 -c1*s2,                                       0;
                0,                                        c2,                                            0;
                -c1*c2,                               s1*s2,                                       0;
                c1*s3 + c3*s1*s2,                -c1*c2*c3,                c3*s1+c1*s2*s3;
                0,                                      -c3*s2,                                 -c2*s3;
                c1*c3*s2-s1*s3,                c2*c3*s1,                c1*c3-s1*s2*s3;
                c1*c3-s1*s2*s3,                c1*c2*s3,                c1*c3*s2-s1*s3;
                0,                                       s2*s3,                              -c2*c3;
                -c3*s1-c1*s2*s3,                -c2*s1*s3,                -c1*s3-c3*s1*s2];
        end;
        
    case 'ZXY',
        R = [c1*c3-s1*s2*s3,                -c2*s1,                c1*s3+c3*s1*s2;
            c3*s1+c1*s2*s3,                c1*c2,                s1*s3-c1*c3*s2;
            -c2*s3,                               s2,                                c2*c3];
        
        if nargout>1,
            dRdvec = [-c3*s1-c1*s2*s3,                -c2*s1*s3,                -c1*s3-c3*s1*s2;
                c1*c3-s1*s2*s3,                c1*c2*s3,                c1*c3*s2-s1*s3;
                0,                                          s2*s3,                             -c2*c3;
                -c1*c2,                               s1*s2,                                  0;
                -c2*s1,                              -c1*s2,                                 0;
                0,                                            c2,                                      0;
                c1*c3*s2-s1*s3,                c2*c3*s1,                c1*c3-s1*s2*s3;
                c1*s3+c3*s1*s2,                -c1*c2*c3,                c3*s1+c1*s2*s3;
                0,                                        -c3*s2,                           -c2*s3];
        end;
        
    case 'ZYX',
        R = [c1*c2,                c1*s2*s3-c3*s1,                s1*s3+c1*c3*s2;
            c2*s1,                c1*c3+s1*s2*s3,                c3*s1*s2-c1*s3;
            -s2,                               c2*s3,                               c2*c3];
        
        if nargout>1,
            dRdvec = [-c2*s1,                    -c1*s2,                                 0;
                c1*c2,                                    -s1*s2,                                0;
                0,                                             -c2,                                    0;
                -c1*c3-s1*s2*s3,                c1*c2*s3,                s1*s3+c1*c3*s2;
                c1*s2*s3-c3*s1,                c2*s1*s3,                c3*s1*s2-c1*s3;
                0,                                          -s2*s3,                                 c2*c3;
                c1*s3-c3*s1*s2,                c1*c2*c3,                c3*s1-c1*s2*s3;
                s1*s3+c1*c3*s2,                c2*c3*s1,                -c1*c3-s1*s2*s3;
                0,                                          -c3*s2,                                -c2*s3];
        end;
        
    otherwise,
        error('Unexpected input for the 2nd argument!');
end;

return;
        

function [vec, dvecdR] = mat2euler(R, op)

bigeps = 1e+5*eps;
[U,~,V] = svd(R);
R = U*V';       % The closest rotation matrix

switch op,
    case 'XYX',
        a2 = acos(R(1,1));
        assert(sin(a2)>bigeps,'The 1st and 3rd Euler angle cannot be determined because of gimbal lock!');
        a1 = atan2(R(2,1),-R(3,1));
        a3 = atan2(R(1,2),R(1,3));
        if nargout>1,
            dvecdR = [0, -R(3,1)/(R(2,1)^2+R(3,1)^2), R(2,1)/(R(2,1)^2+R(3,1)^2),  0, 0, 0,  0, 0, 0;
                -1/(1-R(1,1)^2)^(1/2),    0,    0,    0, 0, 0,     0, 0, 0;
                0,   0,  0, R(1,3)/(R(1,2)^2+R(1,3)^2), 0, 0, -R(1,2)/(R(1,2)^2+R(1,3)^2), 0, 0];
         end;

    case 'XZX',
        a2 = acos(R(1,1));
        assert(sin(a2)>bigeps,'The 1st and 3rd Euler angle cannot be determined because of gimbal lock!');
        a1 = atan2(R(3,1),R(2,1));
        a3 = atan2(R(1,3),-R(1,2));
        if nargout>1,
            dvecdR = [0, -R(3,1)/(R(2,1)^2+R(3,1)^2), R(2,1)/(R(2,1)^2+R(3,1)^2),   0, 0, 0, 0, 0, 0;
                -1/(1-R(1,1)^2)^(1/2),    0,  0,   0, 0, 0,    0, 0, 0;
                0,    0,   0, R(1,3)/(R(1,2)^2 + R(1,3)^2), 0, 0, -R(1,2)/(R(1,2)^2 + R(1,3)^2), 0, 0];
         end;

    case 'YXY',
        a2 = acos(R(2,2));
        assert(sin(a2)>bigeps,'The 1st and 3rd Euler angle cannot be determined because of gimbal lock!');
        a1 = atan2(R(1,2),R(3,2));
        a3 = atan2(R(2,1),-R(2,3));
        if nargout>1,
            dvecdR = [0,  0, 0, R(3,2)/(R(1,2)^2+R(3,2)^2),  0, -R(1,2)/(R(1,2)^2+R(3,2)^2), 0, 0, 0;
                0,  0, 0,   0, -1/(1-R(2,2)^2)^(1/2),  0, 0,   0, 0;
                0, -R(2,3)/(R(2,1)^2+R(2,3)^2), 0,  0,  0,   0, 0, R(2,1)/(R(2,1)^2+R(2,3)^2), 0];
         end;

    case 'YZY',
        a2 = acos(R(2,2));
        assert(sin(a2)>bigeps,'The 1st and 3rd Euler angle cannot be determined because of gimbal lock!');
        a1 = atan2(R(3,2),-R(1,2));
        a3 = atan2(R(2,3),R(2,1));
        if nargout>1,
            dvecdR = [0,  0, 0, R(3,2)/(R(1,2)^2+R(3,2)^2),   0, -R(1,2)/(R(1,2)^2+R(3,2)^2), 0,  0, 0;
                0,    0, 0,    0, -1/(1-R(2,2)^2)^(1/2),   0, 0,    0, 0;
                0, -R(2,3)/(R(2,1)^2+R(2,3)^2), 0,   0,  0,  0, 0, R(2,1)/(R(2,1)^2+R(2,3)^2), 0];
         end;
 
    case 'ZXZ',
        a2 = acos(R(3,3));
        assert(sin(a2)>bigeps,'The 1st and 3rd Euler angle cannot be determined because of gimbal lock!');
        a1 = atan2(R(1,3),-R(2,3));
        a3 = atan2(R(3,1),R(3,2));
        if nargout>1,
            dvecdR = [0, 0,  0, 0, 0,   0, -R(2,3)/(R(1,3)^2+R(2,3)^2), R(1,3)/(R(1,3)^2+R(2,3)^2), 0;
                0, 0,  0, 0, 0,   0,   0,   0, -1/(1-R(3,3)^2)^(1/2);
                0, 0, R(3,2)/(R(3,1)^2+R(3,2)^2), 0, 0, -R(3,1)/(R(3,1)^2+R(3,2)^2), 0,  0,  0];
         end;

    case 'ZYZ',
        a2 = acos(R(3,3));
        assert(sin(a2)>bigeps,'The 1st and 3rd Euler angle cannot be determined because of gimbal lock!');
        a1 = atan2(R(2,3),R(1,3));
        a3 = atan2(R(3,2),-R(3,1));
        if nargout>1,
            dvecdR = [0, 0,  0, 0, 0,  0, -R(2,3)/(R(1,3)^2 + R(2,3)^2), R(1,3)/(R(1,3)^2+R(2,3)^2),  0;
                0, 0,   0, 0, 0,  0,  0,  0, -1/(1-R(3,3)^2)^(1/2);
                0, 0, R(3,2)/(R(3,1)^2+R(3,2)^2), 0, 0, -R(3,1)/(R(3,1)^2 + R(3,2)^2),  0, 0,  0];
         end;

    case 'XYZ',
        a2 = asin(R(1,3));
        assert(cos(a2)>bigeps,'The 1st and 3rd Tait-Bryan angle cannot be determined because of gimbal lock!');
        a1 = atan2(-R(2,3),R(3,3));
        a3 = atan2(-R(1,2),R(1,1));
        if nargout>1,
            dvecdR = [0, 0, 0, 0, 0, 0,  0, -R(3,3)/(R(2,3)^2+R(3,3)^2), R(2,3)/(R(2,3)^2 + R(3,3)^2);
                0, 0, 0,  0, 0, 0, 1/(1 - R(1,3)^2)^(1/2),  0,   0;
                R(1,2)/(R(1,1)^2+R(1,2)^2), 0, 0, -R(1,1)/(R(1,1)^2+R(1,2)^2), 0, 0,  0,   0,   0];
         end;

    case 'XZY',
        a2 = -asin(R(1,2));
        assert(cos(a2)>bigeps,'The 1st and 3rd Tait-Bryan angle cannot be determined because of gimbal lock!');
        a1 = atan2(R(3,2),R(2,2));
        a3 = atan2(R(1,3),R(1,1));
        if nargout>1,
            dvecdR = [0, 0, 0, 0, -R(3,2)/(R(2,2)^2+R(3,2)^2), R(2,2)/(R(2,2)^2+R(3,2)^2),  0, 0, 0;
                0, 0, 0, -1/(1-R(1,2)^2)^(1/2),   0,    0,    0, 0, 0;
                -R(1,3)/(R(1,1)^2+R(1,3)^2), 0, 0,   0,    0,    0, R(1,1)/(R(1,1)^2+R(1,3)^2), 0, 0];
         end;

    case 'YXZ',
        a2 = -asin(R(2,3));
        assert(cos(a2)>bigeps,'The 1st and 3rd Tait-Bryan angle cannot be determined because of gimbal lock!');
        a1 = atan2(R(1,3),R(3,3));
        a3 = atan2(R(2,1),R(2,2));
        if nargout>1,
            dvecdR = [0, 0, 0, 0, 0, 0, R(3,3)/(R(1,3)^2+R(3,3)^2),  0, -R(1,3)/(R(1,3)^2+R(3,3)^2);
                0,  0, 0, 0,  0, 0,  0, -1/(1-R(2,3)^2)^(1/2),   0;
                0, R(2,2)/(R(2,1)^2+R(2,2)^2), 0, 0, -R(2,1)/(R(2,1)^2+R(2,2)^2), 0,  0,  0,  0];
         end;

    case 'YZX',
        a2 = asin(R(2,1));
        assert(cos(a2)>bigeps,'The 1st and 3rd Tait-Bryan angle cannot be determined because of gimbal lock!');
        a1 = atan2(-R(3,1),R(1,1));
        a3 = atan2(-R(2,3),R(2,2));
        if nargout>1,
            dvecdR = [R(3,1)/(R(1,1)^2+R(3,1)^2),   0, -R(1,1)/(R(1,1)^2+R(3,1)^2), 0,  0, 0, 0,  0, 0;
                0, 1/(1-R(2,1)^2)^(1/2),   0, 0,    0, 0, 0,  0, 0;
                0, 0, 0, 0, R(2,3)/(R(2,2)^2+R(2,3)^2), 0, 0, -R(2,2)/(R(2,2)^2+R(2,3)^2), 0];
         end;

    case 'ZXY',
        a2 = asin(R(3,2));
        assert(cos(a2)>bigeps,'The 1st and 3rd Tait-Bryan angle cannot be determined because of gimbal lock!');
        a1 = atan2(-R(1,2),R(2,2));
        a3 = atan2(-R(3,1),R(3,3));
        if nargout>1,
            dvecdR = [0, 0, 0, -R(2,2)/(R(1,2)^2+R(2,2)^2), R(1,2)/(R(1,2)^2+R(2,2)^2),  0, 0, 0, 0;
                0, 0, 0,  0,  0, 1/(1-R(3,2)^2)^(1/2), 0, 0,  0;
                0, 0, -R(3,3)/(R(3,1)^2+R(3,3)^2), 0,  0, 0, 0, 0, R(3,1)/(R(3,1)^2+R(3,3)^2)];
         end;

    case 'ZYX',
        a2 = -asin(R(3,1));
        assert(cos(a2)>bigeps,'The 1st and 3rd Tait-Bryan angle cannot be determined because of gimbal lock!');
        a1 = atan2(R(2,1),R(1,1));
        a3 = atan2(R(3,2),R(3,3));
        if nargout>1,
            dvecdR = [-R(2,1)/(R(1,1)^2+R(2,1)^2), R(1,1)/(R(1,1)^2+R(2,1)^2),  0, 0, 0, 0, 0, 0, 0;
                0,  0, -1/(1-R(3,1)^2)^(1/2), 0, 0,  0, 0, 0,  0;
                0,  0,  0, 0, 0, R(3,3)/(R(3,2)^2+R(3,3)^2), 0, 0, -R(3,2)/(R(3,2)^2+R(3,3)^2)];
         end;

    otherwise,
         error('Unexpected input for the 2nd argument!');
end;

vec = [a1;a2;a3];

return;




%% Test of function:
n = 1000;
om = randn(3,n);
st = {'XYX', 'XZX', 'YXY', 'YZY', 'ZXZ', 'ZYZ', 'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'};
err = zeros(1,n);
for i=1:12,
   for j=1:n,
       R0 = rodrigues(om(:,j));
       R1 = trans_euler_mat(trans_euler_mat(R0, st{i}), st{i});
       err(j) = norm(R1-R0);
   end;
   fprintf(1,['\nError of ' st{i} ' :\n%f\n'],norm(err));
end;

%% Test of dRda:
st = {'XYX', 'XZX', 'YXY', 'YZY', 'ZXZ', 'ZYZ', 'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'};
for j =1:12,
    gain = 1000;
    i=0;
    while gain >100,
        i=i+1;
        om = randn(3,1);
        a = trans_euler_mat(rodrigues(om),st{j});
        da = randn(3,1)/220;
        [R, dRda] = trans_euler_mat(a,st{j});
        R1 = trans_euler_mat(a+da,st{j});
        R2 = R+reshape(dRda*da,3,3);
        gain = norm(R1-R)/norm(R1- R2);   % 大量/小量
    end;
     fprintf(1,['\nFor convention ' st{j} ' :\ni=%d, gain=%f;\n'],[i,gain]);
end;

%% Test of dadR:
st = {'XYX', 'XZX', 'YXY', 'YZY', 'ZXZ', 'ZYZ', 'XYZ', 'XZY', 'YXZ', 'YZX', 'ZXY', 'ZYX'};
for j =1:12,
    gain = 1000;
    i=0;
    while gain >100,
        i=i+1;
        om = randn(3,1);
        R = rodrigues(om);
        dR = rodrigues(om+randn(3,1)/1000)-R;        
        [a, dadR] = trans_euler_mat(R,st{j});
        a1 = trans_euler_mat(R+dR,st{j});
        a2 = a+dadR*dR(:);
        gain = norm(a1-a)/norm(a1- a2);   % 大量/小量
        if gain<1,
            keyboard;
        end;
    end;
     fprintf(1,['\nFor convention ' st{j} ' :\ni=%d, gain=%f;\n'],[i,gain]);
end;