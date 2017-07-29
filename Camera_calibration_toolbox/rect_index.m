function [Irec,ind_new,ind_1,ind_2,ind_3,ind_4,a1,a2,a3,a4] = rect_index(I,R,f,c,k,alpha)
if nargin < 6,
    alpha = 0;
    if nargin < 5,
        k = [0;0;0;0;0];
        if nargin < 4,
            c = [0;0];
            if nargin < 3,
                f = [1;1];
                if nargin < 2,
                    R = eye(3);
                    if nargin < 1,
                        error('ERROR: Need an image to rectify');
                    end;
                end;
            end;
        end;
    end;
end;

KK = [f(1) alpha*f(1) c(1);0 f(2) c(2);0 0 1];

% Note: R is the motion of the points in space
% So: X2 = R*X where X: coord in the old reference frame, X2: coord in the new ref frame.
[nr,nc] = size(I);
Irec = 255*ones(nr,nc);
[mx,my] = meshgrid(1:nc, 1:nr);
npts = nc*nr;
px = reshape(mx,1,npts);
py = reshape(my,1,npts);
xn = KK\[px-1; py-1;ones(1,npts)];

% Rotation: (or affine transformation):
x = R'*xn;
x = [x(1,:)./x(3,:);x(2,:)./x(3,:)];


% Add distortion:
xd = apply_distortion(x,k);
% Reconvert in pixels:
px2 = f(1)*(xd(1,:)+alpha*xd(2,:))+c(1);
py2 = f(2)*xd(2,:)+c(2);

% Interpolate between the closest pixels:
px0 = floor(px2);
py0 = floor(py2);
good_points = find(px0>= 0 & px0<=nx-2 & py0>=0 & py0<=ny-2);

px2 = px2(good_points);
py2 = py2(good_points);
px0 = px0(good_points);
py0 = py0(good_points);

alpha_x = px2 - px0;
alpha_y = py2 - py0;

a1 = (1-alpha_y).*(1-alpha_x);
a2 = (1-alpha_y).*alpha_x;
a3 = alpha_y .* (1-alpha_x);
a4 = alpha_y .* alpha_x;

ind_1 = px0 * nr + py0 + 1;
ind_2 = (px0 + 1) * nr + py0 + 1;
ind_3 = px0 * nr + (py0 + 1) + 1;
ind_4 = (px0 + 1) * nr + (py0 + 1) + 1;

ind_new = (px(good_points)-1)*nr + py(good_points);
Irec(ind_new) = a1 .* I(ind_1) + a2 .* I(ind_2) + a3 .* I(ind_3) + a4 .* I(ind_4);

return;


%% Convert in indices:

fact = 3;
[XX,YY]= meshgrid(1:nc,1:nr);
[XXi,YYi]= meshgrid(1:1/fact:nc,1:1/fact:nr);

%tic;
Iinterp = interp2(XX,YY,I,XXi,YYi); 
%toc

[nri,nci] = size(Iinterp);


ind_col = round(fact*(f(1)*xd(1,:)+c(1)))+1;
ind_row = round(fact*(f(2)*xd(2,:)+c(2)))+1;

good_points = find((ind_col >=1)&(ind_col<=nci)&(ind_row >=1)& (ind_row <=nri));
