function Irect = rectify_image(I,f,c,k,alpha)
% See also distort_image, apply_distortion, comp_distortion.
if nargin < 5,
    alpha = 0;
    if nargin < 4,
        k = [0;0;0;0;0];
        if nargin < 3,
            c = [0;0];
            if nargin < 2,
                f = [1;1];
                if nargin < 1,
                    error('ERROR: Need an image to rectify');
                end;
            end;
        end;
    end;
end;

KK = [f(1), alpha*f(1), c(1);0, f(2), c(2);0, 0, 1];

[ny,nx,nz] = size(I);
npixel = nx*ny;
Irect = I;

% build the index for undistorted pixels
[gx,gy] = meshgrid(1:nx, 1:ny);
px = reshape(gx,1,npixel);       % scan along y direction
py = reshape(gy,1,npixel);
x = KK\[px-1; py-1; ones(1,npixel)];

% Add distortion:
xd = apply_distortion(x(1:2,:),k);
% Reconvert in pixels:
px2 = f(1)*(xd(1,:)+alpha*xd(2,:))+c(1);
py2 = f(2)*xd(2,:)+c(2);

% 2D Interpolation between the four closest pixels:
px0 = floor(px2);       % distorted points
py0 = floor(py2);
good_points = find(px0>= 0 & px0<=nx-2 & py0>=0 & py0<=ny-2);

px2 = px2(good_points);
py2 = py2(good_points);
px0 = px0(good_points);
py0 = py0(good_points);

alpha_x = px2 - px0;
alpha_y = py2 - py0;

clu = (1-alpha_y).*(1-alpha_x);
cru = (1-alpha_y).*alpha_x;
cld = alpha_y .* (1-alpha_x);
crd = alpha_y .* alpha_x;

ind_lu = px0 * ny + py0 + 1;
ind_ru = (px0 + 1) * ny + py0 + 1;
ind_ld = px0 * ny + (py0 + 1) + 1;
ind_rd = (px0 + 1) * ny + (py0 + 1) + 1;

ind_new = (px(good_points)-1)*ny + py(good_points);

for i = 1:nz,
    Iz = I(:,:,i);
    Ir = 255*ones(ny,nx);
    Ir(ind_new) = clu .* Iz(ind_lu) + cru .* Iz(ind_ru) + cld .* Iz(ind_ld) + crd .* Iz(ind_rd);
    Irect(:,:,i) = Ir;
end;

return;


%% Test: inverse function of distort_image
ima_name = './baboon1d.jpg';       %  './D95d1.png';
save_name = './baboon1r.jpg';     % './D95r1.png';
I0 = imread(ima_name);
[nr,nc,nl] = size(I0);
% fov_angle = 35;
% fc = (nc/2)/tan(pi*fov_angle/360) * ones(2,1);
% cc = ([nc;nr]-1)/2;
% div = 2.^(1:5)+5;
% kc = randn(5,1)./div';

% alpha = 0;
I1 = uint8(rectify_image(double(I0),fc,cc,kc));      % ,alpha);
figure(1);
image(I0);
if nl == 1,
    colormap(gray(256));
end;
axis image;
title('Distorted image - Stored in array I0');
figure(2);
image(I1);
if nl == 1,
    colormap(gray(256));
end;
axis image;
title('Rectified image - Stored in array I1');
imwrite(I1,save_name,'png');