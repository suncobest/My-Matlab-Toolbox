% test voxel

vox = outer;
abc = layer(2).abcsize;
nvox = size(vox,2); 
[V,f] = gen_cuboids(vox, abc);

for kk = 1:n_cam
   geom(kk).projverts = project_points_mirror([vox,V], Cam_vec(kk).om, Cam_vec(kk).T, Cam_vec(kk).hand, Cam_vec(kk).fc, Cam_vec(kk).cc, Cam_vec(kk).kc, Cam_vec(kk).alpha_c); 
end

for n_brick = 1:nvox
    figure(1);
    image(im2uint8(bwframe))
    hold on;
    colormap(gray(256));
    set(1,'color',[1 1 1]);
    axis image;
    
    for kk =1:n_cam
        %%% pixelΪһ��ש��8�������ͶӰ��matlab�������� (2,8)
        pixel = geom(kk).projverts(:,nvox+(n_brick-1)*8+1 : nvox+n_brick*8)+1;
        
        %         Pconv = pixel(:,convhull(pixel'));  % ͹���Ķ���������˳�����У��ӵ�һ��������ص���һ����
        %         bwconv = inpolygon(X,Y,Pconv(1,:),Pconv(2,:));  % ����͹���Ķ�ֵͼ
        %         intersection = bwconv & geom(kk).blob;     % ͹����ǰ���Ľ���
        
        % plot one brick in three view
        patch('Faces', f(:,1:6)', 'Vertices', pixel', 'FaceColor', 'b', 'EdgeColor','r','FaceAlpha',0.2);
        plot(geom(kk).projverts(1,n_brick)+1, geom(kk).projverts(2,n_brick)+1, 'y+');  % ����layer(1).voxels�ڸ��ӽǵ�ͶӰ      
    end
    pause;
    hold off;
end
