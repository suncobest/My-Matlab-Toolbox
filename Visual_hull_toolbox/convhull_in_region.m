function [IN, ON] = convhull_in_region( I, vertices, Nvertex, ROI)
% CONVEXHULL_IN_REGION will return the index of convex hull in the
% logical region. IN are logical index of polygon in foreground, while ON
% are logical index of polygon on the border.
%
%   I is the logical image which contain a white foreground region. (ny, nx)
%
%   VERTICES are the vetex points of every convex hull. All convex hulls have
%   the same number of vertices. (2, Np)
%
%   NVERTEX is number of vertices in every convexhull. It's a scalar integer,
%   every convex hull have same number of vertices. So mod(Np, Nvertex)==0.
%   Note: extend the function to line segment, i.e. Nvertex=2
%
%   ROI is the region of interest or boundingbox of the foreground. [x; y; width; height]
%
%   See also inpolygon, convhull, convhulln, convexHull.

assert(islogical(I) && ismatrix(I), 'The 1st input is assumed to be logical image!');
[ny,nx] = size(I);

[m, Np] = size(vertices);
assert(ismatrix(vertices) && m==2,'Unexpected input for the 2nd argument!');

if nargin<3,
    Nvertex = 1;
else
    assert(isscalar(Nvertex) && Nvertex>=1 && Nvertex~=2 && mod(Np,Nvertex)==0, ...
        'Unexpected input for the 3rd argument!');
end;
nhull = Np/Nvertex;
IN = false(1,nhull);
ON = IN;

if nargin==4,
    ROI = ROI(:);
    assert(length(ROI)==4 && all(ROI(1:2)>=0.5) && all(ROI(3:4)>=1) && ROI(1)+ROI(3)<=nx+0.5 ...
        && ROI(2)+ROI(4)<=ny+0.5, 'Unexpected input for the 4th argument!');
    if ~isequal(ROI, [0; 0; nx; ny]+0.5),  % the axis position of image I
        %transform I and vertices
        m = ceil(ROI([2 1]));
        I = I(m(1) : m(1)+ROI(4)-1, m(2) : m(2)+ROI(3)-1);
        vertices = vertices-repmat(ROI(1:2),1,Np)+0.5;
        nx = ROI(3);
        ny = ROI(4);
    end;
end;

if ~any(any(I)),
    return;     % no foreground
end;

if Nvertex==1,
    vertices = round(vertices);
    px = vertices(1,:);
    py = vertices(2,:);
    goodpts = px>=1 & px<=nx & py>=1 & py<=ny;
    ind = py+(px-1)*ny;
    IN(goodpts) = I(ind(goodpts));
    return;
end;

[xxg, yyg] = meshgrid(1:nx, 1:ny);
for kk = 1:nhull,
    ii = (kk-1)*Nvertex;
    px = vertices(1, ii+1:ii+Nvertex)';
    py = vertices(2, ii+1:ii+Nvertex)';
    ind = convhull(px,py);
    in = inpolygon(xxg, yyg, px(ind), py(ind));
    cross = I  & in;
    on = any(any(cross));
    if on,
        in = isequal(cross, in);
        IN(kk) = in;
        ON(kk) = ~in;
    end;
end;

return;



%% Test
flag = 1;   % turn on ROI
Nv = 6;
Npoly = 30;
psize = 10;
nx = 600;
ny = 400;
Im = false(ny,nx);
ind1 = round(ny/5) : round(ny*2/5);
ind2 = round(nx/10) : round(nx/2);
ind3 = round(nx*2/3) : round(nx*5/6);

Im(ind2,ind1)=1;
Im(ind1,ind2)=1;
Im(nx-ind2-round(ny/2), ind3)=1;
Im(ind3-round(ny/3), ind2+round(nx/3))=1;

geom = regionprops(Im);
Bound = geom(1).BoundingBox';
px = nx*rand(1,Npoly);
py = ny*rand(1,Npoly);
px = reshape(px(ones(Nv,1), :),1,[]);
py = reshape(py(ones(Nv,1), :),1,[]);
verts = psize*randn(2,Nv);
verts = repmat(verts, 1, Npoly) +[px;py];

figure(1);
image(Im);
colormap(gray(2));
axis image;
hold on;
plot(verts(1,:), verts(2,:), 'g.');

if flag,
    [inner, border] = convhull_in_region(Im, verts, Nv,Bound);
    rectangle('Position',Bound,'edgecolor','b','linewidth',2);
else
    [inner, border] = convhull_in_region(Im, verts, Nv);
end;
inid = find(inner)
onid = find(border)
if Nv==1,
    plot(verts(1,inid),verts(2,inid),'r+');
    for k = inid,
        xx = verts(:,k);
        text(xx(1),xx(2),num2str(k),'color','m','fontsize',15);
    end;
else
    for k = 1:Npoly,
        i = (k-1)*Nv;
        xx = verts(:,i+1:i+Nv);
        ind = convhull(xx');
        plot(xx(1,ind),xx(2,ind), 'r-');
        cter = mean(xx,2);
        if inner(k),
            text(cter(1),cter(2),num2str(k),'color','m','fontsize',15);
        end;
        if border(k),
            text(cter(1),cter(2),num2str(k),'color','g','fontsize',15);
        end;   
    end;
end;
hold off;

%% texturemap

xg=[0 1;0 1];
yg=[0 0;1 1];
zg=[1 1;1 1];
om=randn(3,1);
T=10*randn(3,1);
x=rodrigues(om)*[xg(:)';yg(:)';zg(:)']+T(:,ones(1,4));
xg(:)=x(1,:);
yg(:)=x(2,:);
zg(:)=x(3,:);
surface(xg, yg, zg,flipdim(uint8(255*repmat(Im,[1,1,3])),1), 'FaceColor','texturemap','EdgeColor','none');

% h= surface('XData',xg,'YData',yg,...
%         'ZData',zg,'CData',flipdim(uint8(255*repmat(Im,[1,1,3])),1),...
%         'FaceColor','texturemap','EdgeColor','none');
