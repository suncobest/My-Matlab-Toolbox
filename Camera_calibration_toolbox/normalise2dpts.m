% NORMALISE2DPTS - normalises 2D homogeneous points
%
% Function translates and normalises a set of 2D homogeneous points
% so that their centroid is at the origin and their mean distance from
% the origin is sqrt(2).  This process typically improves the
% conditioning of any equations used to solve homographies, fundamental
% matrices etc.
%
% Usage:   [newpts, T] = normalise2dpts(pts)
%
% Argument:
%   pts -  3xN array of 2D homogeneous coordinates
%
% Returns:
%   newpts -  3xN array of transformed 2D homogeneous coordinates.  The
%             scaling parameter is normalised to 1 unless the point is at
%             infinity.
%   T      -  The 3x3 transformation matrix, newpts = T*pts
%
% If there are some points at infinity the normalisation transform
% is calculated using just the finite points.  Being a scaling and
% translating transform this will not affect the points at infinity.

% Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/~pk
%
% May 2003      - Original version
% February 2004 - Modified to deal with points at infinity.
% December 2008 - meandist calculation modified to work with Octave 3.0.1
%                 (thanks to Ron Parr)


function [newpts, T] = normalise2dpts(pts)

[m,npts] = size(pts);
if m == 2,
    pts = [pts; ones(1,npts)];
    ind = true(1,npts);
elseif m==3,
    % Find the indices of the points that are not at infinity
    ind = abs(pts(3,:)) > eps;
else
    error('Unexpected dimension of input!');
end


if length(ind) ~= size(pts,2)
    warning('Some points are at infinity');
end

% For the finite points ensure homogeneous coords have scale of 1
pts(1,ind) = pts(1,ind)./pts(3,ind);
pts(2,ind) = pts(2,ind)./pts(3,ind);
pts(3,ind) = 1;

c = mean(pts(1:2,ind),2);      % Centroid of finite points
scale = sqrt(2)/mean(sqrt((pts(1,ind)-c(1)).^2 + (pts(2,ind)-c(2)).^2));

T = [scale   0   -scale*c(1)
      0     scale -scale*c(2)
      0       0      1];

newpts = T*pts;
return;
