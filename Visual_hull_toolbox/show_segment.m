% find the codedir of the toolbox
codeDir = mfilename('fullpath');
flag = find(codeDir=='\' | codeDir=='/',1,'last');
if ~isempty(flag),
    codeDir = codeDir(1 : flag(end)-1);
end;
if exist([codeDir '/load_imgdir.m'],'file'),
    load_imgdir;
    fprintf(1,'\nModify file ''load_imgdir.m'' in code directory if you want to redirect image path!\n');
else
    fprintf(1,'\nNo file ''load_imgdir.m'' found in the code directory!\nPlease compute visual hull first!\n');
    return;
end;

save_name = [imgdir '/segment_variables.mat'];
if exist(save_name,'file')==2,
    fprintf(1,'\nLoading segmentation variable file ''segment_variables.mat'' ...\n');
    load(save_name);
else
    fprintf(1,'\nSegmentation variables ''segment_variables.mat'' not found!\n');
    return;
end;

Nb = 20;
npts = (Nb+1)^2;
[Xe, Ye, Ze] = sphere(Nb);
XYZs = [Xe(:)'; Ye(:)'; Ze(:)'];

ind_active = find(active_images);
N_active = length(ind_active);
assert(isequal(diff(ind_active), ones(1,N_active-1)),'All active images must be continuous to extract kinematics!');
axes_mean = sum(member_axes(:,:,ind_active),3)/N_active;

fprintf(1,'\nThis scipt will show and save segment results of given frames!\n');
if ~exist('show_ellipsoid','var'),
    show_ellipsoid = input('Show the principal ellipsoid for every parts or not? ([]=yes, other=no) ','s');
    show_ellipsoid = isempty(show_ellipsoid);
end;
if ~exist('nlen_vec','var'),
    nlen_vec = input('Length factor of principle axis vector to show: ([]=1.5 times semi_axes) ');
    if isempty(nlen_vec),
        nlen_vec = 1.5;
    end;
end;

while true,
    frame_list = input(['Please input number of frames (1~' num2str(n_frame) '): ']);
    frame_list = round(frame_list);
    if any(frame_list>n_frame | frame_list<1),
        fprintf(1,'\nFrame number out of range (1~%d)!\nPlease input again!\n',n_frame);
        continue;
    end;
    break;
end;

for kk = frame_list,
    if ~active_images(kk),
        warning('Frame %d is not active!',kk);
        continue;
    end;
    count = ceil(kk/nfpu);
    base = (count-1)*nfpu;
    bk = kk-base;
    nc = sprintf(['%0' ndigit 'd'],count);
    save_name = [imgdir '/segment_part' nc '.mat'];
    if exist(save_name,'file')==2,
        load(save_name);
    else
        fprintf(1,'\nERROR: cannot find data ''%s''!\n',save_name);
        return;
    end;
    XX = X_cell{bk};
    seg_id = seg_cell{bk};
    
    % plot segment result
    plot_string = 'plot3(';
    for k=1:nparts+1,
        plot_string = [plot_string 'XX(1,seg_id==' num2str(k) '),XX(3,seg_id==' ...
            num2str(k) '),-XX(2,seg_id==' num2str(k) '),''.'','];
    end;
    plot_string = [plot_string 'XX(1,seg_id==0),XX(3,seg_id==0),-XX(2,seg_id==0),''.'')'];
    figure(3); hold off;
    eval(plot_string); hold on;
    for ii = 1:nparts+1,
        % plot the three semi axes
        ct = member_center(:,ii,kk);
        ax = axes_mean(:,ii);
        vc = rodrigues(member_omega(:,ii,kk));
        ctt = ct(:,ones(1,3));
        xyz = ctt+vc*diag(ax)*nlen_vec;
        for i=1:3,
            plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',2);
        end;
        if show_ellipsoid,
            % generate ellipsoid from unit sphere
            XYZe = vc*diag(ax)*XYZs+ct(:,ones(1,npts));
            Xe(:) = XYZe(1,:);
            Ye(:) = XYZe(2,:);
            Ze(:) = XYZe(3,:);
            surf(Xe,Ze,-Ye,'FaceColor',palette(ii,:), 'EdgeColor',palette(ii,:), 'FaceAlpha',0.05);
        end;
    end;
    axis equal vis3d off;
    view(az,el);
    set(3,'color',[1 1 1]*0.7);
    drawnow;
    framenb = sprintf(['%0' ndigit 'd'],kk);
    print(3,'-djpeg',[imgdir '/segment_' framenb '.jpg']);
    print(3,'-deps',[imgdir '/segment_' framenb '.eps']);
    saveas(3, [imgdir '/segment_' framenb],'fig');
end;

fprintf(1,'\nDone!\n');
