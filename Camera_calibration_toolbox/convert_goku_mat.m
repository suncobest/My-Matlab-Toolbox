% convert
np1D=3;
n_cam=9;
n_ima=498;
imsize=[1280;1024]*ones(1,n_cam);
rodlen = [0, 125, 375];

n_view = n_ima * n_cam;
npts = np1D*n_ima;
hand_list = ones(1,n_cam);
xx = zeros(n_ima,7,n_cam);
for i=1:n_cam,
    xx(:,:,i)=load(sprintf('camera%02d.txt',i-1));
end;
xx(:,1,:)=[];
active_imgviews =  permute(any(xx,2),[3,1,2]);

xx = reshape(permute(xx,[2,3,1]),[2,3,n_view]);
x_cell = cell(1,n_view);
ind_active_views = find(active_imgviews(:)');
for kth = ind_active_views,
    x_cell{kth} = xx(:,:,kth);
end;
active_images = sum(active_imgviews,1)>1;
active_imgviews(:,~active_images) = 0;
ind_active = find(active_images);

xx=reshape(xx,[2,3,n_cam,n_ima]);
xx=reshape(permute(xx,[1,2,4,3]),[2,npts,n_cam]);

string_save = ['save multicam_real_data active_imgviews active_images npts n_ima n_cam ' ...
                   'n_view np1D rodlen hand_list imsize xx x_cell'];
eval(string_save);
