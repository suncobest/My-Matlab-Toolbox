function [ out ] = quatlog( in )
% QUATLOG - The natural logarithm of quaternions.  In are quaternions with
%   4*n elements. 1~3 rows are imaginary part (3d vector part), the 4th row
%   are real part (scalor part). If IN = [x, y, z, w], T=norm([x, y, z, w]);
%   [x, y, z]/T = stheta*[x, y, z]/norm([x, y, z]), ctheta=w/T;
%   then out = [theta*[x, y, z]/norm([x, y, z]), log(T)].
%
%   Bibliography: Quaternions, Interpolation and Animation, p15~16.
%   See also quatexp, quatmul.

%  by zpf, form BIT, 2015-7-22

if isvector(in),
    in = in(:);
end;
assert(ismatrix(in) && size(in,1)==4,'The argument must be a matrix with 4 rows!');

in2 = in.*in;
T = sqrt(sum(in2));
out = in./T(ones(4,1),:);
theta = acos(out(4,:)); 
in2 = sqrt(sum(in2(1:3,:)));
id = in2~=0;
theta(id) = theta(id)./in2(id);
out = [in(1:3,:).*theta(ones(3,1),:); log(T)];
% turning theta around n=[x,y,z] is the same rotation as turning 2*pi-theta
% around the anti-axis -n. As a result, the 1st quaternion is the same as
% the 2nd one up to a minus sign.

return;


%% Test
p = sym('p%d',[4,1]);
q = sym('q%d',[4,1]);
p = sym(p,'real');
q = sym(q,'real');
r = quatlog([p,q]);
r = simplify(r)

%%
a = randn(4,10);
err = quatexp(quatlog(a))-a
err = quatlog(quatexp(a))-a
