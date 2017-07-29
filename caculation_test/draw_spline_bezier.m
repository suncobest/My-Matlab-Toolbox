%% Demo: spline_bezier
filename = 'spline_bezier.gif';
lw1 = 3.0;
lw2 = 1.5;

Np = 5;
ndim = 2;
n2 = 80;

Nd1 = ones(ndim,1);
x=1:Np;
xx=linspace(1, Np, n2);
y =10*rand(ndim,Np);
yy = spline_bezier(x,y,xx);

hi = diff(x);
% caculate the auxiliary points q
Q = y;
a= hi(2:Np-1)./(2*(hi(1:Np-2)+hi(2:Np-1)));
b = 0.5-a;
Q(:, 2:Np-1) = y(:,2:Np-1)*3/2-a(Nd1,:).*y(:,1:Np-2)-b(Nd1,:).*y(:,3:Np);

a = [y yy Q];
min_yy = min(a,[],2);
max_yy = max(a,[],2);
delta = max(abs(a(:)))/60;
interval  = [min_yy(1)-0.5, max_yy(1)+0.5, min_yy(2)-0.5, max_yy(2)+0.5];

figure(1);
if ndim==2,
    plot(y(1,:),y(2,:),'co', y(1,:),y(2,:),'c-', Q(1,:),Q(2,:),'g+', Q(1,:),Q(2,:),'g-', yy(1,:),yy(2,:),'b-','linewidth', lw2);
elseif ndim==3,
    plot3(y(1,:),y(2,:),y(3,:), 'co', y(1,:),y(2,:),y(3,:), 'c-', Q(1,:), Q(2,:), Q(3,:),'g+', ...
        Q(1,:), Q(2,:), Q(3,:), 'g-', yy(1,:),yy(2,:),yy(3,:),'b-', 'linewidth', lw2);
end;

%%
rgb = lines(4);
dt2 = 0.2;
el = 20;
[~,idx] = histc(xx,[-inf,x(2:Np-1),inf]);
tt = (xx-x(idx))./hi(idx);

for t=1:n2,
    % draw
    figure(1);
    hold off;
    if ndim==2,
        plot(y(1,:),y(2,:),'c-', Q(1,:),Q(2,:),'g-',yy(1,:),yy(2,:),'b-', 'linewidth', lw1);
    elseif ndim==3,
        plot3(y(1,:),y(2,:),y(3,:),'c-', Q(1,:), Q(2,:), Q(3,:), 'g-', yy(1,:),yy(2,:),yy(3,:),'b-', 'linewidth', lw1);
    end;
    hold on;
    for i=1:Np,
        if ndim==2,
            plot(y(1,i),y(2,i),'c+',Q(1,i),Q(2,i),'g+', 'markersize',10.0,'linewidth', lw2);
            text(y(1,i)+delta,y(2,i)+delta,['p', num2str(i)],'fontsize',16,'fontweight','bold');
            text(Q(1,i)+delta,Q(2,i)+delta,['q', num2str(i)],'fontsize',16,'fontweight','bold');
        elseif ndim==3,
            plot3(y(1,i),y(2,i),y(3,i),'c+',Q(1,i), Q(2,i), Q(3,i),'g+', 'markersize',10.0,'linewidth', lw2);
            text(y(1,i)+delta,y(2,i)+delta,y(3,i)+delta,['p', num2str(i)],'fontsize',16,'fontweight','bold');
            text(Q(1,i)+delta,Q(2,i)+delta,Q(3,i)+delta,['q', num2str(i)],'fontsize',16,'fontweight','bold');
        end;
    end;
    
    ti = tt(t);
    ii = idx(t);
    Ai = (1-ti)*y(:, ii)+ti*y(:,ii+1);
    Bi = (1-ti)*Q(:, ii)+ti*Q(:,ii+1);
    ti = 2*ti*(1-ti);
    yi = (1-ti)*Ai+ti*Bi;
    if norm(yi-yy(:,t))>1e-10,
        error('Result not match!');
    end;
    if ndim==2,
        plot([Ai(1);Bi(1)], [Ai(2);Bi(2)],'color', rgb(4,:),'linewidth', lw2);
        plot(Ai(1),Ai(2),'mo', Bi(1),Bi(2),'ro', yi(1),yi(2),'ko','markersize',6.0,'linewidth', lw2);
        axis(interval);
    elseif ndim==3,
        plot3([Ai(1);Bi(1)], [Ai(2);Bi(2)],[Ai(3);Bi(3)], 'color', rgb(4,:),'linewidth', lw2);
        plot3(Ai(1),Ai(2), Ai(3),'mo', Bi(1),Bi(2),Bi(3),'ro', yi(1),yi(2),yi(3),'ko','markersize',6.0,'linewidth', lw2);
        axis([interval, min_yy(3)-0.5, max_yy(3)+0.5]);
        view(t*360/n2, el);
    end;
    
    set(gcf,'color','w');
    set(gcf, 'unit', 'normalized', 'position', [0.1,0.1,0.8,0.8]);  % 设置窗口全屏
    set(gca,'position',[0 0 1 1]); % 设置绘图区域（坐标系）在窗口中的位置, position = [left, bottom, width, height]
    axis off equal;
    
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [I,map] = rgb2ind(im,256);
    
    if  t == 1,
        imwrite(I,map,filename,'gif','LoopCount',Inf,'DelayTime',dt2);
    else
        imwrite(I,map,filename,'gif','WriteMode','append','DelayTime',dt2);
    end;
end;