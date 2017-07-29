function [GX, GY] = gaussGradient(Im, sigma)

% Create an even-length 1-D separable Derivative of Gaussian filter

% Determine filter length
filterLength = 8*ceil(sigma);
n = (filterLength - 1)/2;
x = -n:n;

% Create 1-D Gaussian Kernel
ssq=sigma*sigma;
gau = exp(-x.^2/(2*ssq));     % 一维高斯滤波器

% Normalize to ensure gauss filter sums to one
sumgau = sum(gau);
if sumgau~=0,
    gau = gau/sumgau;
end;

% Create 1-D Derivative of Gaussian Kernel
dgaudx = -x.*gau/ssq;  %高斯滤波器的一阶导数
negVals = dgaudx < 0;
posVals = dgaudx > 0;
dgaudx(posVals) = dgaudx(posVals)/sum(dgaudx(posVals));
dgaudx(negVals) = dgaudx(negVals)/abs(sum(dgaudx(negVals)));

% Compute smoothed numerical gradient of image I along x (horizontal)
% direction. GX corresponds to dG/dx, where G is the Gaussian Smoothed
% version of image I.
GX = imfilter(Im, gau', 'conv','replicate');
GX = imfilter(GX, dgaudx, 'conv','replicate');

% Compute smoothed numerical gradient of image I along y (vertical)
% direction. GY corresponds to dG/dy, where G is the Gaussian Smoothed
% version of image I.
GY = imfilter(Im, gau, 'conv','replicate');
GY = imfilter(GY, dgaudx', 'conv','replicate');
