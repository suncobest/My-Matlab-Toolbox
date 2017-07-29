function h = plot_lines_std(x,y,deltay,varargin)
% FUNCTION plot_lines_std plot lines along with standard deviation.
[x,y,deltay,style,palette,lw,xt,yt,sw] = parse_inputs(x,y,deltay,varargin);
[m,n] = size(x);

fs = 9;      % font size of latex mode
fn = 'Arial';   % font name of latex mode
d = 0.1;    % margin of figure
w = 1-1.5*d;    % width and height of figure
trans = 0.3;    % alpha to set fill transparency
width = 5;      % unit: inch
height = 4;
resolution = 600;   % dpi
save_name = 'mean_std';

% xlim and ylim
a = x(:);
a = round([min(a),max(a)]);
xx = [x; x(m:-1:1,:)];
yy = [y+deltay; y(m:-1:1,:)-deltay(m:-1:1,:)];
b = [min(yy(:)), max(yy(:))];
b = round(b+[-1,1]*(b(2)-b(1))/10);

h = figure('renderer','opengl'); 
axes('position',[d d w w]);
hold on;
for i=1:n,
    plot(x(:,i), y(:,i), style{i}, 'color', palette(:,i)', 'linewidth',lw(i));
    fill(xx(:,i),yy(:,i),palette(:,i)','EdgeColor','none','FaceAlpha',trans);
end;

if ~isempty(xt),
    if sw(1),
        xlabel(xt, 'Interpreter','latex','fontname',fn,'fontsize',fs);
    else
        xlabel(xt);
    end;
end;
if ~isempty(yt),
    if sw(2),
        ylabel(yt, 'Interpreter','latex','fontname',fn,'fontsize',fs);
    else
        ylabel(yt);
    end;
end;

set(gca,'XMinorTick','on','YMinorTick','on','box','off','xlim',a,'ylim',b);
set(h,'color','w','PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 width height]);
drawnow;
print(h,'-djpeg',['-r' num2str(resolution)],save_name);
saveas(gcf, save_name,'fig');

return;


%-------------------------------------------------------
function [x,y,deltay,style,palette,lw,xt,yt,sw] = parse_inputs(x,y,deltay,v)
%PARSE_INPUTS
%   [x,y,deltay,style,palette,lw,xt,yt,sw] = parse_inputs(x,y,deltay,v) returns line style, colors,
%   width along with xlabel text xt, ylabel text yt, latex interpreter switch sw.
%   x: m*n;
%   y: m*n;
%   deltay:, m*n;
%   style: cell of 1*n; like: {'-', ':', '-.', '--'}
%   palette: numeric color map of 3*n, or char color vector of 1*n;
%   lw: width of lines, 1*n;
%   xt: xlabel text;
%   yt: ylabel text;
%   sw: latex switch of xlabel and ylabel;

[m1,n1] = size(y);
[m,n] = size(x);
if m==1,
    m = n;
    n = n1;
    x = x(:)*ones(1,n);
end;
if n==1,
    n = n1;
    x = x(:,ones(1,n));
end;
assert(isequal([m,n],[m1,n1]) && isequal(size(deltay),[m,n]),['Dimension of inputs do not match ' ...
                    'with each other!']);

style  = cell(1,n);
palette = lines(n)';
for i=1:n,
    style{i} = '-';
end;
lw = ones(1,n);
xt = [];
yt = [];
sw = zeros(1,2);

if ~isempty(v),
    nv = length(v);
    n1 = length(v{1});      % n1 is the length of cell vector
    if ischar(v{1}),            % style like '--', ':', or '-.' is an element, not cell any more
        n1 = 1;
    end;
    switch n1,
        case 0,
            % do nothing
        case 1,
            for i =1:n,
                style{i} = v{1};
            end;
        case n,
            style = v{1};
        otherwise,
            error('Number of styles do not match with lines!');
    end;
    
    if nv>1,
        [m1,n1] = size(v{2});
        if n1,
            flag = 0;
            if ischar(v{2}),
                flag = 1;
                assert(m1==1,'Unexpected input for colors!');
            else
                assert(m1==3, 'The color matrix must have 3 rows!');
            end;
            if n1==1,
                palette = v{2}*ones(1,n);
            elseif n1==n,
                palette = v{2};
            else
                error('Number of colors do not match with lines!');
            end;
            if flag,
                palette = char(palette);
            end;
        end;
        if nv>2,
            n1 = length(v{3});
            switch n1,
                case 0,
                    % do nothing
                case 1,
                    lw = v{3}*ones(1,n);
                case n,
                    lw = v{3};
                otherwise,
                    error('Number of line specification do not match with lines!');
            end;
            if nv>3,
                if ~isempty(v{4}),
                    xt = v{4};
                end;
                if nv>4,
                    if ~isempty(v{5}),
                        yt = v{5};
                    end;
                    if nv>5 && ~isempty(v{6}),
                        n1 = length(v{6});
                        switch n1,
                            case 0,
                                % do nothing
                            case 1,
                                sw = v{6}*ones(1,2);
                            case 2,
                                sw = v{6};
                            otherwise,
                                error('Unexpected dimension of latex switch!');
                        end;
                    end;
                end;
            end;
        end;
    end;
end;

return;



%% Test
m = 4;
n = 10;
x = (0:n)/n;
y = randn(m,n+1);
dy = randn(m,n+1)/10;
x1 = linspace(0,1,200);
y1 = spline_interp3(x,y,x1);
dy1 = spline_interp3(x,dy,x1);
% plot_lines_std(x1,y1',dy1',{':','-','-.','--'},'rgbm', [0.2, 0.5, 0.7, 1],...
%     '\it\fontname{Arial}\fontsize{9}x','\it\fontname{Arial}\fontsize{9}y')
plot_lines_std(x1,y1',dy1',[],[], 0.2,'$\hat t$', '$\theta^*\ \ (^\circ)$',[1 1])
