%%  spline_interp3-bezier curve gif 
%  clamp the curves with 0 tangent (jointing bezier splines)
filename = 'Jointing_bezier_curve.gif';
lw1 = 3.0;
lw2 = 1.5;

Np = 8;
rgb = lines(Np);

ndim = 2;
x=1:Np;
n2 = 80;
xx=linspace(1, Np, n2);
y =10*rand(ndim,Np);

[yy, coef] = spline_interp3(x,y,xx,'c');    % bezier curve-clamped spline with zero boundary tangents.
hi = diff(x);
ci = coef(:, 3);
ci = reshape(ci,ndim,[]);
Ai = ci/3+y(:,1:Np-1);
Bi = y(:,2:Np);
ci = hi(1:Np-2)./hi(2:Np-1);
Bi(:,1:Np-2) = Bi(:,1:Np-2) +(Bi(:,1:Np-2) -Ai(:, 2:Np-1)).*ci(ones(ndim,1),:);
ci = [y(:,1:Np-1); Ai; Bi];
yab = [reshape(ci,ndim,[]),y(:,Np)];
n=size(yab,2);
min_yy = min(yab,[],2);
max_yy = max(yab,[],2);
delta = max(abs(yab(:)))/50;
interval  = [min_yy(1)-0.5, max_yy(1)+0.5, min_yy(2)-0.5, max_yy(2)+0.5];
figure(1);
if ndim==2,
    plot(y(1,:),y(2,:),'ro', Ai(1,:),Ai(2,:),'y+', Bi(1,:),Bi(2,:),'y+', yy(1,:),yy(2,:),'b-','linewidth', lw2);
elseif ndim==3,
    plot3(y(1,:),y(2,:),y(3,:),'ro', Ai(1,:),Ai(2,:),Ai(3,:),'y+', Bi(1,:),Bi(2,:),Bi(3,:),'y+', yy(1,:),yy(2,:),yy(3,:),'b-','linewidth', lw2);
end;

%% subdivide bezeier
hold on;
y2=subdivision_bezier(yab,1);
if ndim ==2,
    plot(y2(1,:), y2(2,:),'mo:', yab(1,:), yab(2,:),'y:');
    for i = 1:n,
        text(yab(1,i)+delta, yab(2,i)+delta,  num2str(i), 'color', 'r');
    end;
elseif ndim==3,
    plot3(y2(1,:), y2(2,:), y2(3,:),'mo:', yab(1,:), yab(2,:),yab(3,:), 'y:');
    for i = 1:n,
        text(yab(1,i)+delta, yab(2,i)+delta, yab(3,i)+delta, num2str(i), 'color', 'r');
    end;
end;


%%
dt2 = 0.2;
el = 20;
[~,idx] = histc(xx,[-inf,x(2:Np-1),inf]);
tt = (xx-x(idx))./hi(idx);
for t=1:n2,
    % draw
    figure(1);
    hold off;
    if ndim==2,
        plot(yab(1,:),yab(2,:),'b-','linewidth', lw1);  % 'color',rgb(1,:)
    elseif ndim==3,
        plot3(yab(1,:),yab(2,:),yab(3,:),'b-','linewidth', lw1);  % 'color',rgb(1,:)
    end;
    hold on;
    for i=1:n,
        if ndim==2,
            plot(yab(1,i),yab(2,i),'r+','markersize',10.0,'linewidth', lw2);
            text(yab(1,i)+delta,yab(2,i)+delta,num2str(i),'fontsize',16,'fontweight','bold');
        elseif ndim==3,
            plot3(yab(1,i),yab(2,i),yab(3,i),'r+','markersize',10.0,'linewidth', lw2);
            text(yab(1,i)+delta,yab(2,i)+delta,yab(3,i)+delta,num2str(i),'fontsize',16,'fontweight','bold');
        end;
    end;
    if ndim==2,
        plot(yy(1,:),yy(2,:),'g-','linewidth', lw1);  % 'color',rgb(2,:)
    elseif ndim==3,
        plot3(yy(1,:),yy(2,:),yy(3,:),'g-','linewidth', lw1);  % 'color',rgb(2,:)
    end; 
    
    ti = tt(t);
    ii = idx(t);
    yi = [y(:, ii), Ai(:,ii),Bi(:,ii),y(:,ii+1)]; 
    for k = 4:-1:2,                                      % De Casteljau algorithm
        yi = (1-ti)*yi(:,1:k-1)+ti*yi(:,2:k);
        if ndim==2,
            plot(yi(1,:),yi(2,:),'color',rgb(k,:),'linewidth', lw2);
        elseif ndim==3,
            plot3(yi(1,:),yi(2,:),yi(3,:),'color',rgb(k,:),'linewidth', lw2);
        end;
    end;
    if ndim==2,
        plot(yi(1),yi(2),'ko','markersize',5.0,'linewidth', lw2);
        axis(interval);
    elseif ndim==3,
        plot3(yi(1),yi(2),yi(3),'ko','markersize',5.0,'linewidth', lw2);
        axis([interval, min_yy(3)-0.5, max_yy(3)+0.5]);
        view(t*360/n2, el);
    end;
        
    set(gcf,'color','w');
    set(gcf, 'unit', 'normalized', 'position', [0.1,0.1,0.8,0.8]);  % 设置窗口全屏
    set(gca,'position',[0 0 1 1]); % 设置绘图区域（坐标系）在窗口中的位置, position = [left, bottom, width, height]
    axis off equal;
    
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [I,map] = rgb2ind(im,256);
    
    if  t == 1,
		imwrite(I,map,filename,'gif','LoopCount',Inf,'DelayTime',dt2);
	else
		imwrite(I,map,filename,'gif','WriteMode','append','DelayTime',dt2);
    end;  
end;