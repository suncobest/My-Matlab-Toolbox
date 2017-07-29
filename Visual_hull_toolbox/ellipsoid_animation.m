% ellipsoid animation
h = figure;
for kk=1:nn,
    figure(h);hold off;
    for ii =1:nparts+1,
        ax = axes_mean(:,ii);
        ct = ctposition(:,ii,kk);
        vc = rodrigues(axisAngle(:,ii,kk));
        % plot axes
        ctt = ct(:,ones(1,2));
        xyz = ctt+vc(:,1:2)*diag(ax(1:2))*nlen_vec;
        for i=1:2,
            plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',2);
            hold on;
        end;
        % generate ellipsoid from unit sphere
        XYZe = vc*diag(ax)*XYZs+ct(:,ones(1,npts));
        patch('Faces', faces, 'Vertices', [XYZe(1,:); XYZe(3,:); -XYZe(2,:)]','FaceColor',palette(ii,:),...
            'EdgeColor',palette(ii,:),'FaceAlpha',0.1);
    end;
    axis equal; grid on;
    set(gcf, 'renderer', 'zbuffer','color',[1,1,1]*0.7);
    cameratoolbar('ResetCameraAndSceneLight');
    cameratoolbar('Togglescenelight');
    view(az,el);
    axis(axes_range);
    pause(0.001);
end;
