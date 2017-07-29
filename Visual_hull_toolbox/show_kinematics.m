%% show orientation kinematics
nx = 1280;
ny = 800;
height = 4;

a = [ts(1),ts(end)]*1e3;
d = 0.1;    % margin of figure
w = 1-d*1.5;    % width and height of figure
w1 = 0.5;
w2 = 0.5;
e = a(1)-(a(2)-a(1))*d/(2*w);   % x postion of ylabel
% color: 1 blue; 2 green; 3 red; 4 cyan; 5 magenta; 6 yellow; 7 black;
% if oldflag==1, old data (blue,green,red), reprojected new data (cyan,magenta,yellow)
palette = lines(7);
m = 12;     % font size

fprintf(1,'\nShow insect kinematics (position and orientation) in subplots.\n');
for n=4:5,
    h = (w-(n-1)*d/2)/n;
    c = h+d/2;
    c = (n-4:n-1)*c+d;
    figure;
    if n==5,
        % draw body center position
        if oldflag,
            XX = zeros(6,nn);
            XX(1:3,:) = reshape(ctposition1(:,1,:),3,[]);
            XX(1:3,:) = R1'*(XX(1:3,:)-XX(1:3,1)*ones(1,nn));
            XX(4:6,:) = reshape(ctposition(:,1,:),3,[]);
            XX(4:6,:) = R0'*(XX(4:6,:)-XX(4:6,1)*ones(1,nn));
            b = [min(XX(:)), max(XX(:))];
            delta = (b(2)-b(1))/10;
            b = round(b+[-1,1]*delta);
            axes('position',[d d w h]);
            plot(1e3*ts, XX(4,:), '-', 'color', palette(5,:), 'linewidth', w1);   % magenta
            hold on;
            plot(1e3*ts, XX(5,:), '-', 'color', palette(4,:), 'linewidth', w1);   % cyan
            plot(1e3*ts, XX(6,:), '-', 'color', palette(6,:), 'linewidth', w1);   % yellow
            plot(1e3*ts, XX(1,:), '-', 'color', palette(1,:), 'linewidth', w2);   % blue
            plot(1e3*ts, XX(2,:), '-', 'color', palette(2,:), 'linewidth', w2);   % green
            plot(1e3*ts, XX(3,:), '-', 'color', palette(3,:), 'linewidth', w2);   % red
        else
            XX = reshape(ctposition(:,1,:),3,[]);
            XX = R0'*(XX-XX(:,1)*ones(1,nn));
            b = [min(XX(:)), max(XX(:))];
            delta = (b(2)-b(1))/10;
            b = round(b+[-1,1]*delta);
            axes('position',[d d w h]);
            plot(1e3*ts, XX,'linewidth',w1);
        end;
        text(1e3*ts(end),XX(1,end)+delta,['\it\fontname{Arial}\fontsize{' num2str(m) '}\color[rgb]{', ...
            num2str(palette(1,:)), '}X_b'],'HorizontalAlignment','center');
        text(1e3*ts(end),XX(2,end)+delta,['\it\fontname{Arial}\fontsize{' num2str(m) '}\color[rgb]{', ...
            num2str(palette(2,:)), '}Y_b'],'HorizontalAlignment','center');
        text(1e3*ts(end),XX(3,end)+delta,['\it\fontname{Arial}\fontsize{' num2str(m) '}\color[rgb]{', ...
            num2str(palette(3,:)), '}Z_b'],'HorizontalAlignment','center');
        set(gca,'YMinorTick','on','box','off','xlim',a,'ylim',b);
        xlabel(['\fontname{Arial}\fontsize{' num2str(m) '}time (ms)']);
        ylabel(['\fontname{Arial}\fontsize{' num2str(m) '}center (mm)'], 'position',[e, sum(b)/2]);
        %     legend('Xb','Yb','Zb');
    end;
    
    % draw body orientation kinematics
    if oldflag,
        XX = [yaw1; pitch1; roll1; yaw; pitch; roll];
    else
        XX = [yaw; pitch; roll; yaw1; pitch1; roll1];
    end;
    b = [min(XX(:)), max(XX(:))];
    delta = (b(2)-b(1))/10;
    b = round(b+[-1,1]*delta);
    axes('position',[d c(1) w h]);
    plot(1e3*ts, XX(4,:), '-', 'color', palette(5,:), 'linewidth', w1);   % magenta
    hold on;
    plot(1e3*ts, XX(5,:), '-', 'color', palette(4,:), 'linewidth', w1);   % cyan
    plot(1e3*ts, XX(6,:), '-', 'color', palette(6,:), 'linewidth', w1);   % yellow
    plot(1e3*ts, XX(1,:), '-', 'color', palette(1,:), 'linewidth', w2);   % blue
    plot(1e3*ts, XX(2,:), '-', 'color', palette(2,:), 'linewidth', w2);   % green
    plot(1e3*ts, XX(3,:), '-', 'color', palette(3,:), 'linewidth', w2);   % red
    text(1e3*ts(end),XX(1,end)+delta,['\fontname{Arial}\fontsize{' num2str(m) '}\color[rgb]{', ...
        num2str(palette(1,:)), '}Yaw'],'HorizontalAlignment','center');
    text(1e3*ts(end),XX(2,end)+delta,['\fontname{Arial}\fontsize{' num2str(m) '}\color[rgb]{', ...
        num2str(palette(2,:)), '}Pitch'],'HorizontalAlignment','center');
    text(1e3*ts(end),XX(3,end)+delta,['\fontname{Arial}\fontsize{' num2str(m) '}\color[rgb]{', ...
        num2str(palette(3,:)), '}Roll'],'HorizontalAlignment','center');
    if n==5,
        set(gca,'XTick',[],'xcolor','w','YMinorTick','on','box','off','xlim',a,'ylim',b);
    else
        set(gca,'YMinorTick','on','box','off','xlim',a,'ylim',b);
        xlabel(['\fontname{Arial}\fontsize{' num2str(m) '}time (ms)']);
    end;
    ylabel(['\fontname{Arial}\fontsize{' num2str(m) '}angle (^\circ)'], 'position',[e, sum(b)/2]);
    % legend('yaw','pitch','roll');
    
    % draw wing orientation kinematics
    if oldflag,
        XX = [strpsi1; strpsi];
    else
        XX = [strpsi; strpsi1];
    end;
    b = [min(XX(:)), max(XX(:))];
    delta = (b(2)-b(1))/10;
    b = round(b+[-1,1]*delta);
    axes('position',[d c(2) w h]);
    plot(1e3*ts, XX(3,:), '-', 'color', palette(4,:), 'linewidth', w1);   % right wing: cyan
    hold on;
    plot(1e3*ts, XX(4,:), '-', 'color', palette(6,:), 'linewidth', w1);   % left wing: yellow
    plot(1e3*ts, XX(1,:), '-', 'color', palette(2,:), 'linewidth', w2);   % right wing: green
    plot(1e3*ts, XX(2,:), '-', 'color', palette(3,:), 'linewidth', w2);   % left wing: red
    ylabel(['\fontname{Arial}\fontsize{' num2str(m) '}incidence (^\circ)'],'position',[e, sum(b)/2]);
    set(gca,'XTick',[],'xcolor','w','YMinorTick','on','box','off','xlim',a,'ylim',b);
    % legend('R','L');
    
    if oldflag,
        XX = [strtheta1; strtheta];
      else
        XX = [strtheta; strtheta1];
    end;
    b = [min(XX(:)), max(XX(:))];
    delta = (b(2)-b(1))/10;
    b = round(b+[-1,1]*delta);
    axes('position',[d c(3) w h]);
    plot(1e3*ts, XX(3,:), '-', 'color', palette(4,:), 'linewidth', w1);   % right wing: cyan
    hold on;
    plot(1e3*ts, XX(4,:), '-', 'color', palette(6,:), 'linewidth', w1);   % left wing: yellow
    plot(1e3*ts, XX(1,:), '-', 'color', palette(2,:), 'linewidth', w2);   % right wing: green
    plot(1e3*ts, XX(2,:), '-', 'color', palette(3,:), 'linewidth', w2);   % left wing: red
    ylabel(['\fontname{Arial}\fontsize{' num2str(m) '}deviation (^\circ)'], 'position',[e, sum(b)/2]);
    set(gca,'XTick',[],'xcolor','w','YMinorTick','on','box','off','xlim',a,'ylim',b);
    % legend('R','L');
    
    if oldflag,
        XX = [strphi1; strphi];
    else
        XX = [strphi; strphi1];
    end;
    b = [min(XX(:)), max(XX(:))];
    delta = (b(2)-b(1))/10;
    b = round(b+[-1,1]*delta);
    axes('position',[d c(4) w h]);
    plot(1e3*ts, XX(3,:), '-', 'color', palette(4,:), 'linewidth', w1);   % right wing: cyan
    hold on;
    plot(1e3*ts, XX(4,:), '-', 'color', palette(6,:), 'linewidth', w1);   % left wing: yellow
    plot(1e3*ts, XX(1,:), '-', 'color', palette(2,:), 'linewidth', w2);   % right wing: green
    plot(1e3*ts, XX(2,:), '-', 'color', palette(3,:), 'linewidth', w2);   % left wing: red
    ylabel(['\fontname{Arial}\fontsize{' num2str(m) '}amplitude (^\circ)'],'position',[e, sum(b)/2]);
    set(gca,'XTick',[],'xcolor','w','YMinorTick','on','box','off','xlim',a,'ylim',b);
    % legend('R','L');
    
    % time of reversal
    if ind_up(1) < ind_down(1),
        ind = [ind_up(1:length(ind_down)); ind_down];
    else
        ind = [ind_up(1:length(ind_down)-1); ind_down(2:end)];
    end;
    if ind_up(end)>ind(end),
        ind = [ind, [ind_up(end); ind_active(end)]];
    end;
    b = (ind-ind_active(1))/FramePS*1e3;	% period of upstoke
    % plot transparent color bar on upstrokes
    axes('position',[d, d, w, w]); hold on;
    for ii=1:size(b,2),
        temp = [b(1,ii), b(2,ii), b(2,ii), b(1,ii); 0, 0, 1, 1];
        fill(temp(1,:),temp(2,:),[1 1 1]*0.2,'EdgeColor','none','FaceAlpha',0.3);
    end;
    set(gca,'xlim',a); axis off;
    resolution = round(ny/height);   % dpi
    set(gcf,'color','w','PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]*2/resolution);
    drawnow;
    save_name = [imgdir '/kinematics' num2str(n-4)];
    print(gcf,'-djpeg',['-r' num2str(resolution*2)],save_name);
    saveas(gcf, save_name,'fig');
end;
