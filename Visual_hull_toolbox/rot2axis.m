function [ x,y,z ] = rot2axis( az,el )
% ROT2AXIS(X) change the 'view direction' format from a composition angle 
% '[az,el]' to a 3d vector 'x'. Vector 'x' is the view axis, normalized 
% as a unit vector here.The compsition '[az,el]' is short for azimuth and 
% elevation respectively.
% Take azimuth and elevation as the first two Euler angles 'z-x',which 
% denotes two successive rotations around the Euler axes 'z' and 'x''.
% Looking through the 'y'' axis is The view direction. 
%
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
% [x,y,z] is a unit vector, so its magnitude equals 1.
% 要求角度全部为degree度

if nargin==2
    az=[az,el];
end

if ~isvector(az) || length(az)~=2 || ~isa(az,'double')  % limit az
    error('Input must be a 2D vector !');
else
    x = cosd(az(2)) * sind(az(1));
    y = -cosd(az(2)) * cosd(az(1));
    z = sind(az(2));
end

switch nargout  % 3 output arguments
    case 0
        x=[x,y,z];
    case 1
        x=[x,y,z];
    case 2
        error('Wrong numbers of output arguments!');
end

end

