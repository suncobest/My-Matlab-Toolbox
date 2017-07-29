function xd = add_distortion(xu,k_dist,center,fc)
% Add distortion to image.
% k_dist is the distortion factor
% center is the center of distortion, fc is the focal length

if nargin < 4,
    fc = [1;1];
    if nargin < 3,
        center = [0;0];
        if nargin < 2 || norm(k_dist)==0,
            xd = xu;
            return;
        end;
    end;
end;
[m,n] = size(xu);
assert(m==2,'Unexpected points dimension!');
assert(length(k_dist)<6 && ~isempty(k_dist),'Unexpected dimension of distortion coefficient!');

if length(fc)==1,
    fc = fc*[1;1];
end;

center = center(1:2);

x_centered = (xu-center*ones(1,n))./(fc*ones(1,n));
r2 = x_centered(1,:).^2 + x_centered(2,:).^2;

k = zeros(5,1);
k(k_dist~=0) = k_dist(k_dist~=0);
k_dist = k;
if norm(k(2:5))==0,
    k_dist = k(1);
    xd = center*ones(1,n) + (x_centered.*(1+k_dist*ones(2,1)*r2)).*(fc*ones(1,n));
    return;
end;


r4 = r2.^2;
r6 = r2.^3;
% Radial distortion:
rdist = 1 + k_dist(1) * r2 + k_dist(2) * r4 + k_dist(5) * r6; 

xd = x_centered .* (ones(2,1)*rdist);

% tangential distortion:
a1 = 2*x_centered(1,:).*x_centered(2,:);
a2 = r2 + 2*x_centered(1,:).^2;
a3 = r2 + 2*x_centered(2,:).^2;

delta_x = [k_dist(3)*a1 + k_dist(4)*a2 ;
    k_dist(3) * a3 + k_dist(4)*a1];

xd = xd + delta_x;

xd = xd.*(fc*ones(1,n))+ center*ones(1,n);

return;



%% test
Nx = 10;
Ny = 8;
edge = 0.2;
fc = max([Nx-1;Ny-1]*edge)/3;
[Ygrid,Xgrid]=meshgrid(0:edge:(Ny-1)*edge, 0:edge:(Nx-1)*edge);
center = edge*[Nx-1;Ny-1]/2;   % edge*[Nx-1;Ny-1]/2;
x0=[Xgrid(:)';Ygrid(:)'];
k = [-0.3;0.1;0;0.05;]; % -0.3;
xdist = add_distortion(x0,k,center,fc);
xgrid_d = zeros(Nx,Ny);
ygrid_d = zeros(Nx,Ny);
xgrid_d(:) = xdist(1,:);
ygrid_d(:) = xdist(2,:);

figure(1);
plot(x0(1,:),x0(2,:),'r.')            % original data
hold on,axis equal
plot(Xgrid,Ygrid,'r-',Xgrid',Ygrid','r-')

plot(xdist(1,:),xdist(2,:),'g.')            % distorted data
plot(xgrid_d,ygrid_d,'g-',xgrid_d',ygrid_d','g-');

set(gcf,'color','w');
set(gca,'ydir','reverse');
axis tight


