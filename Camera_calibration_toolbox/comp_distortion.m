function x_comp  = comp_distortion(x_dist,k2)

%       [x_comp] = comp_distortion(x_dist,k2);
%
%       compensates the radial distortion of the camera
%       on the image plane.
%
%       x_dist : the image points got with radial distortion.
%       x_comp : The image plane points after correction for the distortion
%
%       x_comp and x_dist are 2xN arrays
%
%       NOTE : This compensation has to be done after the substraction
%              of the center of projection, and division by the focal
%              length.
%
%       (do it up to a third order approximation)


assert(size(x_dist,1)==2,'The dimension of the points should be 2xN');
if length(k2) > 1,
    x_comp  = comp_distortion_oulu(x_dist,k2);
else
    x_comp = x_dist;
    for i = 1:3,
        radius_2= sum(x_comp.^2);
        radial_distortion = 1 + k2 * ones(2,1)* radius_2;
        x_comp = x_dist./radial_distortion;
    end;
end;

%% Function completely checked : It works fine !!!
