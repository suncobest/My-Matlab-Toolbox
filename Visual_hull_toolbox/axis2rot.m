function [ az, el ] = axis2rot( x, y, z)
% AXIS2ROT(X) change the view-direction format from a 3d vector 'x' to a 
% composition of angles '[az,el]', short for azimuth and elevation respectively.
%
% Take azimuth and elevation as the first two Euler angles 'z-x',which 
% denotes two successive rotations around the Euler axes 'z' and 'x''.
% Looking through the 'y'' axis is the view direction. 

% Both the vector 'x' and the two angles '[az,el]'can be used as
% arguments for the function 'VIEW' to set the view direction.
% For example:
% view(x)
% view(az,el)
%
%   Detailed explanation goes here
% -y axis correspond to the composition angles [0,0]. 
% Define 'Axis <> Rot' as a pair of identical 'view direction'. 
% Then we have: 
% [0,-1,0] <> [0,0];
% [0,1,0] <> [180,0];
% [0,0,1] <> [0,90]; 
% [0,0,-1] <> [0,-90];
% [1,0,0] <> [90,0];
% [-1,0,0] <> [-90,0];
%
% 角度全部化为degree度

switch nargin % 3 arguments
    case 3
        x=[x,y,z];
    case 2
        x=[x,y];    
end

if ~isvector(x) || length(x)~=3 || norm(x)==0  % limit x
    error('Input must be a nonzero 3D vector !');
else
   
%     if size(x,1)==3  % check if x is a column vector
%         x3=sqrt(x'*x);
%         x2=sqrt(x(1:2)'*x(1:2));
%     else  % or x is a row vector
%         x3=sqrt(x*x');
%         x2=sqrt(x(1:2)*x(1:2)');
%     end 
    x3=norm(x); 
    x2=norm(x(1:2));
    el=asind(x(3)/x3);
    if x2==0
        az=0;
    else
        if x(2)>=0  % acosd among [0,180]
            az=90+acosd(x(1)/x2);
            if az>180
                az=az-360;  % az among (-180,180]
            end
        else
            az=90-acosd(x(1)/x2);
        end
    end
end

if nargout<2
    az=[az,el];
end

end

