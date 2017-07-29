function [xc,good,outspeed] = cornertracker(xt,I0,I1,wintx,winty,wxtrack,wytrack,inspeed)

% [xc] = cornertracker(xt,I);
%
% Finds the sub-pixel corners on the image I with initial guess xt (in the previous image)
% xt and xc are 2xN matrices. The first component is the x coordinate
% (horizontal) and the second component is the y coordinate (vertical)
%
% first convert to image matrix index i and j.
% wintx is radius of window in the horizontal direction; (j)
% winty is radius of window in the vertical direction; (i)
%
% wintx and winty are the radius of neighbourhood to extract conners.
% wxtrack and wytrack are the radius of neighbourhood to track conners.
%
% Based on Lucas and kanade's algorithm and Harris corner finder method (cornerfinder).
%
% Finds corners to a precision below .1 pixel!
% good is a binary vector indicating wether a feature point has been properly found.
%
% see cornerfinder.m

% Beijing Institute of Technology
% zpf -- Jan. 27th, 2015

if nargin < 5,
    wintx = 7;
    winty = 7;
end;

if nargin < 8,
    inspeed = zeros(size(xt));
    if nargin < 7,
        wxtrack = 3*wintx;
        wytrack = 3*winty;
    end;
end;

xt = fliplr(xt'); % 取x方向为行标i增大的方向, y方向为列标j增加的方向，则画面法向朝外
inspeed = fliplr(inspeed');

% the grid after transformation: x-->j, y-->i
mask = exp(-((-wintx:wintx)'/(wintx)).^2) * exp(-((-winty:winty)/(winty)).^2);

dx_tolerance = 0.5;  % motion tolerance
resolution = 0.005;  % corner points extraction
MaxIter = 20;

offx = [-wintx:wintx]'*ones(1,2*winty+1);
offy = ones(2*wintx+1,1)*[-winty:winty];

% shiftx = [-wxtrack:wxtrack]'*ones(1,2*wytrack+1);
% shifty = ones(2*wxtrack+1,1)*[-wytrack:wytrack];
% N = size(xt,1);

[nx,ny] = size(I0);
good = all(xt ~= 0,2); % 1: good features points for tacking; 0: lost
idx = find(good)';
xt = xt.*[good,good];
inspeed = inspeed.*[good,good];
outspeed = inspeed;

xc0 = xt + inspeed; % first guess...
xc = xc0;
for i = idx,
    cIx = xt(i,1); 			%
    cIy = xt(i,2); 			% Coords. of the point
    crIx = round(cIx); 		% on the initial image
    crIy = round(cIy); 		%
    itIx = cIx - crIx; 		% Coefficients to compute the sub pixel accuracy.
    itIy = cIy - crIy;
    if itIx > 0,            % cIx > crIx,  cIx:(0, 0.5)
        vIx = [itIx 1-itIx 0]'; 	% mask for interpolation
    else                    % cIx <= crIx,  cIx:[0.5, 0]
        vIx = [0 1+itIx -itIx]';
    end;
    if itIy > 0,
        vIy = [itIy 1-itIy 0];
    else
        vIy = [0 1+itIy -itIy];
    end;

    % What if the sub image is not in?
    if (crIx-wxtrack-2 < 1),
        xmin=1; xmax = 2*wxtrack+5;
    elseif (crIx+wxtrack+2 > nx),
        xmax = nx; xmin = nx-2*wxtrack-4;
    else
        xmin = crIx-wxtrack-2; xmax = crIx+wxtrack+2;
    end;

    if (crIy-wytrack-2 < 1),
        ymin=1; ymax = 2*wytrack+5;
    elseif (crIy+wytrack+2 > ny),
        ymax = ny; ymin = ny-2*wytrack-4;
    else
        ymin = crIy-wytrack-2; ymax = crIy+wytrack+2;
    end;

    SI0 = double(I0(xmin:xmax,ymin:ymax)); % The necessary neighborhood
    SI0 = conv2(conv2(SI0,vIx,'same'),vIy,'same');  % The subpixel interpolated neighborhood (对窗口内各点进行线性插值，得到网格平移到位置[cIx,cIy]处的像素)
    SI0 = SI0(2:2*wxtrack+4,2:2*wytrack+4); % omit the boundary

%     SI1 = double(I1(xmin:xmax,ymin:ymax));
%     SI1 = conv2(conv2(SI1,vIx,'same'),vIy,'same');
%     SI1 = SI1(2:2*wxtrack+4,2:2*wytrack+4);
%     gt = SI1-SI0;                         % the difference of I0 and I1
%     gt = gt(2:2*wxtrack+2,2:2*wytrack+2); % omit the boundary

    [gy,gx] = gradient(SI0); 		% The gradient image
    gx = gx(2:2*wxtrack+2,2:2*wytrack+2); % extraction of the useful parts of the gradients
    gy = gy(2:2*wxtrack+2,2:2*wytrack+2); % omit the boundary

    gxx = gx .* gx;
    gyy = gy .* gy;
    gxy = gx .* gy;

%     bb = -[sum(sum(gx .* gt)); sum(sum(gy .* gt))];
    a = sum(sum(gxx));
    b = sum(sum(gxy));
    c = sum(sum(gyy));
    G = [a b;b c];    % G*dx=bb

    dx = dx_tolerance + 1;  % just larger than dx_tolerance

    % corner tracking iteration
    compt = 0; 				% no iteration yet
    while (norm(dx) > dx_tolerance) && (compt<MaxIter),
        cIx = xc(i,1);
        cIy = xc(i,2); 			% Coords. of the point
        crIx = round(cIx); 		% on the initial image
        crIy = round(cIy);
        itIx = cIx - crIx; 		% Coefficients to compute the sub pixel accuracy.
        itIy = cIy - crIy;
        if itIx > 0,            % cIx > crIx,  cIx:(0, 0.5)
            vIx = [itIx 1-itIx 0]'; 	% mask for interpolation
        else                    % cIx <= crIx,  cIx:[0.5, 0]
            vIx = [0 1+itIx -itIx]';
        end;
        if itIy > 0,
            vIy = [itIy 1-itIy 0];
        else
            vIy = [0 1+itIy -itIy];
        end;

        % What if the sub image is not in?
        if (crIx-wxtrack-2 < 1),
            xmin=1; xmax = 2*wxtrack+5;
        elseif (crIx+wxtrack+2 > nx),
            xmax = nx; xmin = nx-2*wxtrack-4;
        else
            xmin = crIx-wxtrack-2; xmax = crIx+wxtrack+2;
        end;

        if (crIy-wytrack-2 < 1),
            ymin=1; ymax = 2*wytrack+5;
        elseif (crIy+wytrack+2 > ny),
            ymax = ny; ymin = ny-2*wytrack-4;
        else
            ymin = crIy-wytrack-2; ymax = crIy+wytrack+2;
        end;

        SI1 = double(I1(xmin:xmax,ymin:ymax));
        SI1 = conv2(conv2(SI1,vIx,'same'),vIy,'same');
        SI1 = SI1(2:2*wxtrack+4,2:2*wytrack+4);
        gt = SI1-SI0;                         % the difference of I0 and I1
        gt = gt(2:2*wxtrack+2,2:2*wytrack+2); % omit the boundary

        bb = -[sum(sum(gx .* gt)); sum(sum(gy .* gt))];
        dx = (G\bb)';   % motion of feature point
        xc(i,:) = [cIx,cIy] + dx;
        compt = compt + 1;
    end;

    xc2 = xc(i,:);
    dx = xc2-xc0(i,:);
    if abs(dx(1))>wxtrack || abs(dx(2))>wytrack || xc2(1)<1 || xc2(1)> nx || xc2(2)<1 || xc2(2)>ny,
%         keyboard;
        xc(i,:) = [0 0];
        good(i) = 0;
        outspeed(i,:) = [0 0];
        continue;
    end;

    % corner extracting iteration
    v_extra = resolution + 1; 		% just larger than resolution
    compt = 0; 				% no iteration yet
    while (norm(v_extra) > resolution) && (compt<MaxIter),
        cIx = xc(i,1); 			
        cIy = xc(i,2); 			% Coords. of the point
        crIx = round(cIx); 		% on the initial image
        crIy = round(cIy); 		
        itIx = cIx - crIx; 		% Coefficients to compute the sub pixel accuracy.
        itIy = cIy - crIy;
        if itIx > 0,
            vIx = [itIx 1-itIx 0]'; 	% mask for interpolation
        else
            vIx = [0 1+itIx -itIx]';
        end;
        if itIy > 0,
            vIy = [itIy 1-itIy 0];
        else
            vIy = [0 1+itIy -itIy];
        end;

        % What if the sub image is not in?
        if (crIx-wintx-2 < 1),
            xmin=1; xmax = 2*wintx+5;
        elseif (crIx+wintx+2 > nx),
            xmax = nx; xmin = nx-2*wintx-4;
        else
            xmin = crIx-wintx-2; xmax = crIx+wintx+2;
        end;

        if (crIy-winty-2 < 1),
            ymin=1; ymax = 2*winty+5;
        elseif (crIy+winty+2 > ny),
            ymax = ny; ymin = ny-2*winty-4;
        else
            ymin = crIy-winty-2; ymax = crIy+winty+2;
        end;

        SI1 = double(I1(xmin:xmax,ymin:ymax)); % The necessary neighborhood
        SI1 = conv2(conv2(SI1,vIx,'same'),vIy,'same');  % The subpixel interpolated neighborhood (对窗口内各点进行线性插值，得到网格平移到位置[cIx,cIy]处的像素)
        SI1 = SI1(2:2*wintx+4,2:2*winty+4); % omit the boundary
        [gy,gx] = gradient(SI1); 		% The gradient image
        gx = gx(2:2*wintx+2,2:2*winty+2); % extraction of the useful parts of the gradients
        gy = gy(2:2*wintx+2,2:2*winty+2); % omit the boundary

        px = cIx + offx;   % 网格平移到位置[cIx,cIy]处的新网格坐标
        py = cIy + offy;

        gxx = gx .* gx .* mask;
        gyy = gy .* gy .* mask;
        gxy = gx .* gy .* mask;

        % 平均的角点坐标: 利用自相关矩阵[gxx, gxy; gxy, gyy]对位置进行加权平均
        % bb = sum([gxx, gxy; gxy, gyy]*[px; py]) = sum([gxx, gxy; gxy, gyy])* [X; Y]
        % 其中xc2 =[X; Y], G  = sum([gxx, gxy; gxy, gyy]) = [a, b;b, c]，则有
        % xc2 = inv(G) * bb = adjacent(G)*bb/det(G) = [c, -b;-b,a]*bb/(a*c-b^2)

        bb = [sum(sum(gxx .* px + gxy .* py)); sum(sum(gxy .* px + gyy .* py))];

        a = sum(sum(gxx));
        b = sum(sum(gxy));
        c = sum(sum(gyy));
        dt = a*c - b^2;
        xc2 = [c*bb(1)-b*bb(2) a*bb(2)-b*bb(1)]/dt;

        v_extra = xc(i,:) - xc2;   % 当两次迭代的像素差长度norm(v_extra)<=resolution时，认为已收敛，跳出循环
        xc(i,:) = xc2;
        compt = compt + 1;
    end;

    dx = xc2-xc0(i,:);
    if abs(dx(1))>wxtrack || abs(dx(2))>wytrack,
        xc(i,:) = [0 0];
        good(i) = 0;
        outspeed(i,:) = [0 0];
        continue;
    end;
    outspeed(i,:) = xc2-xt(i,:);
end;

xc = fliplr(xc);
xc = xc';
outspeed = fliplr(outspeed);
outspeed = outspeed';
good = good';

