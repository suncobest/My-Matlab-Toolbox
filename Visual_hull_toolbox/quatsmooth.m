function [ Qm ] = quatsmooth( Qs, Ws, SW )
% QUATSMOOTH - smooth a quaternion sequnence like image filter. Qs is a quaternion sequnence
%   with 4*n elements like n pixels. 1~3 rows are imaginary parts (3d vector part), the 4th row
%   is real part (scalor part). Ws is the weighted filter to perform row wise 'correlation' on
%   the input quaternion sequnence. SW is the output format switch to turn on 'same' scale matrix.
%   If SW==1, then Qm is of the same size as input Qs, else Qm is the valid part of smooth result.
%   If length(Ws) is even, then SW=0.
%
%   reference: Markley, Averaging Quaternions, Journal of guidance, control, and
%   dynamics, 2007, p1194.

%  by zpf, form BIT, 2016-12-25

[m,n] = size(Qs);
assert(ismatrix(Qs) && m==4,'The 1st argument must be a matrix with 4 rows!');
assert(isvector(Ws) && all(Ws>=0),'The 2nd argument must be a non-negative vector');
lw = length(Ws);
if nargin<3,
    SW = true;
else
    SW = ~~SW;
end;
if ~mod(lw,2),  % l is even, interpolate in the middle
    SW = false;
end;

Qs2 = Qs.*Qs;
T = sqrt(sum(Qs2));
Qs = Qs./T(ones(4,1),:);
m = n-lw+1;
if SW,
    Qm = Qs;
    kk = (lw-1)/2;
else
    n = floor((lw-1)/2); 
    Qm = Qs(:,n+1:n+m);
    kk = 0;
end;

Qs2 = zeros(4,4,lw);
for ii=1:m,
    for jj=1:lw,
        T = Qs(:,ii+jj-1);
        Qs2(:,:,jj) = Ws(jj)*(T*T');
    end;
    temp = sum(Qs2,3);
    [~,~,T] = svd(temp);
    if sum(T(:,1).*Qm(:,kk+ii))<0,
        Qm(:,kk+ii) = -T(:,1);
    else
        Qm(:,kk+ii) = T(:,1);
    end;
end;

return;


%% Test
n = 20;
delta = 1/20;
x0= 1:n;
xx=1:.05:n;
lx = length(xx);
y=randn(3,n);
y = y./(ones(3,1)*sqrt(sum(y.*y)));
p1 = y(:,1);
for i=2:n,
    q1 = y(:, i);
    if sum(p1.*q1)<0,
        p1= -q1;
        y(:, i) = p1;
    else
        p1 = q1;
    end;
end;
y0 = [zeros(1,n); y];
yy0 = squad(x0,y0,xx);
l = 9;
f = ones(1,l)/l;
% f = fspecial('gaussian',[1,l],2);
yyn = yy0+[zeros(1,lx);randn(3,lx)/20];
yyn = yyn./(ones(4,1)*sqrt(sum(yyn.*yyn)));
yys = quatsmooth(yyn,f);

hf = 1;
figure(hf);
set(gcf,'color','w');
h = plot3(y0(2,:), y0(3,:), y0(4,:), 'go', yy0(2,:), yy0(3,:), yy0(4,:), 'g-', ...
    yyn(2,:), yyn(3,:), yyn(4,:),'m.', yys(2,:), yys(3,:), yys(4,:), 'r-');
axis equal off;
axis([-1  1 -1 1 -1 1]);
% zdir = [0 0 1];
% rotate(h,zdir,90)     % 旋转句柄对象

for i = 1:n,
    text(y0(2,i)+delta, y0(3,i)+delta, y0(4,i)+delta, num2str(i), 'color', 'g');
end;

% filename = ['quatsmooth' num2str(hf)];
filename = 'quatsmooth';
view(3);
% [az,el]=view;
az0 = 0;
el = 20;
dt = 0.1;  % gif 的帧时间
% 沿az方向旋转一周
daz = 5;  % 每一帧方位角（azimuth）的增量
n = 360/daz;   % 旋转一周步数

for i = 1 : n-1,
    view(az0 + daz*(i-1),el);
    %     rotate(h,zdir,daz) ;    % 旋转句柄对象
    drawnow;
    frame = getframe(hf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    if  i == 1,
        imwrite(A,map,[filename,'.gif'],'gif','LoopCount',Inf,'DelayTime',dt);
    else
        imwrite(A,map,[filename,'.gif'],'gif','WriteMode','append','DelayTime',dt);
    end;
end;
% saveas(gcf, filename, 'fig');
