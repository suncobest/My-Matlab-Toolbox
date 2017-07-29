function [ out ] = quatexp( in )
% QUATEXP - The exponentiation of quaternions.  In are quaternions with 4*n
%   elements. 1~3 rows are imaginary part (3d vector part), the 4th row are
%   real part (scalor part). If IN = [x, y, z, w], theta=norm([x, y, z]);
%   then out = exp(w)*[[x, y, z]*sin(theta)/theta, cos(theta)].
%
%   Bibliography: Quaternions, Interpolation and Animation, p15~16.
%   See also quatlog, quatmul.

%  by zpf, form BIT, 2015-7-22

bigeps = 1e-7;

if isvector(in),
    in = in(:);
end;
assert(ismatrix(in) && size(in,1)==4,'The argument must be a matrix with 4 rows!');

vec = in(1:3,:);
theta = sqrt(sum(vec.*vec));
index = theta>=bigeps;
out = theta(index);
if any(index),
    stheta= sin(out)./out;                      % sin(eps)/eps =1
    vec(:,index) = vec(:,index).*stheta(ones(3,1),:);
end;
out = [vec; cos(theta)];
theta = exp(in(4,:));
out = theta(ones(4,1),:).*out;

return;


%% Test
a = randn(4,10);
err = quatexp(quatlog(a))-a
err = quatlog(quatexp(a))-a

%%
p = sym('p%d',[4,1]);
q = sym('q%d',[4,1]);
p = sym(p,'real');
q = sym(q,'real');
r = quatexp([p,q]);
r = simplify(r)

