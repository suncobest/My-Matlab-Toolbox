% COMPUTE_FUNDMATRIX - computes fundamental matrix from 8 or more points
%
% Function computes the fundamental matrix from 8 or more matching points in
% a stereo pair of images.  The normalised 8 point algorithm given by
% Hartley and Zisserman p265 is used.  To achieve accurate results it is
% recommended that 12 or more points are used
%
% Usage:   [F, e1, e2] = compute_fundmatrix(x1, x2)
%          [F, e12] = compute_fundmatrix(x)
%
% Arguments:
%          x1, x2 - Two sets of corresponding 3xN set of homogeneous
%          points.
%
%          x      - If a single argument is supplied it is assumed that it
%                   is in the form x = [x1; x2]
% Returns:
%          F      - The 3x3 fundamental matrix such that x2'*F*x1 = 0.
%          e1     - The epipole in image 1 such that F*e1 = 0
%          e2     - The epipole in image 2 such that F'*e2 = 0
%

% Copyright (c) 2002-2005 Peter Kovesi


function [F,e1,e2] = compute_fundmatrix(varargin)
[x1, x2, npts] = checkargs(varargin(:));

% Normalise each set of points so that the origin
% is at centroid and mean distance from origin is sqrt(2).
% normalise2dpts also ensures the scale parameter is 1.

[x1, T1] = normalise2dpts(x1);
[x2, T2] = normalise2dpts(x2);

% Build the constraint matrix
A = [x2(1,:)'.*x1(1,:)',  x2(1,:)'.*x1(2,:)',  x2(1,:)', ...
     x2(2,:)'.*x1(1,:)',  x2(2,:)'.*x1(2,:)',  x2(2,:)', ...
     x1(1,:)',             x1(2,:)',           ones(npts,1)];

[U,D,V] = svd(A,0);

% Extract fundamental matrix from the column of V corresponding to
% smallest singular value.
F = reshape(V(:,9),3,3)';

% Enforce constraint that fundamental matrix has rank 2 by performing
% a svd and then reconstructing with the two largest singular values.
[U,D,V] = svd(F,0);
F = U*diag([D(1,1) D(2,2) 0])*V';

% Denormalise
F = T2'*F*T1;

if nargout > 1,  	% Solve for epipoles
    [U,D,V] = svd(F,0);
    if V(3,3)>eps,
        e1 = V(:,3)/V(3,3);
    end;
    if U(3,3)>eps,
        e2 = U(:,3)/U(3,3);
    end;
    if nargout==2,
        e1 = [e1,e2];
    end;
end;

return;


%--------------------------------------------------------------------------
% Function to check argument values and set defaults
function [x1, x2, npts] = checkargs(arg);

if length(arg) == 2,
    x1 = arg{1};
    x2 = arg{2};
    if ~all(size(x1)==size(x2)),
        error('x1 and x2 must have the same size');
    elseif size(x1,1) ~= 3,
        error('x1 and x2 must be 3xN');
    end;

elseif length(arg) == 1,
    if size(arg{1},1) ~= 6,
        error('Single argument x must be 6xN');
    else
        x1 = arg{1}(1:3,:);
        x2 = arg{1}(4:6,:);
    end;
else
    error('Wrong number of arguments supplied');
end;

npts = size(x1,2);
if npts < 8,
    error('At least 8 points are needed to compute the fundamental matrix');
end;

return;
