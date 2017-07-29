filename = 'bezier_curve.gif';

Np = 8;
ndim = 2;
yy=10*rand(Np,ndim);

%%
min_yy = min(yy);
max_yy = max(yy);
delta = max(abs(yy(:)))/50;

n=size(yy,1);
ndiv = 30;
b=bezier_interp(yy',ndiv)';    % function is transposed to have every point in the column direction

rgb = lines(Np);
dt2 = 0.2;
dt = 1/ndiv;

el = 20;

for t = 0:dt:1,
    
    % draw
    figure(1);
    hold off;
    if ndim==2,
        plot(yy(:,1),yy(:,2),'b-','linewidth',3.0);  % 'color',rgb(1,:)
    elseif ndim==3,
        plot3(yy(:,1),yy(:,2),yy(:,3),'b-','linewidth',3.0);  % 'color',rgb(1,:)
    end;
    hold on;
    for i=1:n,
        if ndim==2,
            plot(yy(i,1),yy(i,2),'r+','markersize',10.0,'linewidth',1.5);
            text(yy(i,1)+delta,yy(i,2)+delta,num2str(i),'fontsize',16,'fontweight','bold');
        elseif ndim==3,
            plot3(yy(i,1),yy(i,2),yy(i,3),'r+','markersize',10.0,'linewidth',1.5);
            text(yy(i,1)+delta,yy(i,2)+delta,yy(i,3)+delta,num2str(i),'fontsize',16,'fontweight','bold');
        end;
    end;
    if ndim==2,
        plot(b(:,1),b(:,2),'g-','linewidth',3.0);  % 'color',rgb(2,:)
    elseif ndim==3,
        plot3(b(:,1),b(:,2),b(:,3),'g-','linewidth',3.0);  % 'color',rgb(2,:)
    end;
    
    yi = yy;
    for k = Np:-1:2,                                      % De Casteljau algorithm
        yi = (1-t)*yi(1:k-1,:)+t*yi(2:k,:);
        if ndim==2,
            plot(yi(:,1),yi(:,2),'color',rgb(k,:),'linewidth',1.5);
        elseif ndim==3,
            plot3(yi(:,1),yi(:,2),yi(:,3),'color',rgb(k,:),'linewidth',1.5);
        end;
    end;
    if ndim==2,
        plot(yi(1),yi(2),'ko','markersize',5.0,'linewidth',1.5);
        axis([min_yy(1)-0.5, max_yy(1)+0.5, min_yy(2)-0.5, max_yy(2)+0.5]);
    elseif ndim==3,
        plot3(yi(1),yi(2),yi(3),'ko','markersize',5.0,'linewidth',1.5);
        axis([min_yy(1)-0.5, max_yy(1)+0.5, min_yy(2)-0.5, max_yy(2)+0.5,min_yy(3)-0.5, max_yy(3)+0.5]);
        view(t*360,el);
    end;
        
    set(gcf,'color','w');
    set(gcf, 'unit', 'normalized', 'position', [0.1,0.1,0.8,0.8]);  % 设置窗口全屏
    set(gca,'position',[0 0 1 1]); % 设置绘图区域（坐标系）在窗口中的位置, position = [left, bottom, width, height]
    axis off equal;
    
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [I,map] = rgb2ind(im,256);
    
    if  t == 0,
		imwrite(I,map,filename,'gif','LoopCount',Inf,'DelayTime',dt2);
	else
		imwrite(I,map,filename,'gif','WriteMode','append','DelayTime',dt2);
    end;  
end;
