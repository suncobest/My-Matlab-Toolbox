function draw_arc(AX,center,radius,dstart,dend,ndiv)
% draw an arc at the center of the radius from the dstart to the dend
% ndiv is the arc subdivision number, AX is the handle number of figure;
if nargin <6,
    ndiv = 100;
    if nargin < 5,
        dstart =0;
        dend = 360;
        if nargin<3,
            radius=1;
            if nargin<2,
                center=[0,0];
                if nargin<1,
                    AX=1;
                end;
            end;
        end;
    end;
end;

if ~(isvector(center) && length(center)==2),
    error('The second argument (arc center) must be a 2D vector!');
end;

if ~isscalar(radius),
    error('The third argument (arc radius) must be a scalar!');
end;

if ~isscalar(dstart),
    error('The fourth argument (arc start angle) must be a scalar!');
end;

if ~isscalar(dend),
    error('The fifth argument (arc end angle) must be a scalar!');
end;

if ~isscalar(ndiv),
    error('The sixth argument (arc subdivision number) must be a scalar!');
end;

theta = linspace(dstart,dend,ndiv);
x = center(1)+radius*cosd(theta);
y = center(2)+radius*sind(theta);
figure(AX);
axis image;
plot(x,y,'b','linewidth',2);

return;