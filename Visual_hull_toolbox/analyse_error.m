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

save_name = [imgdir '/visualhull_environment.mat'];
if exist(save_name,'file')==2,
    fprintf(1,'\nLoading common data from file ''visualhull_environment.mat'' ...\n');
    load(save_name);
else
    fprintf(1,'\n''visualhull_environment.mat'' not found! You must run script ''init_visualhull'' first!\n');
    return;
end;

if  false,
    n = 5;
    MaxIter = 10;
    n3 = n^3;
    center3d = subvoxcenter(center3d_mat, n, boundsize);
    center2d = NaN(2, n3*n_frame, n_cam);
    if calib_mode,
        for pp = 1:n_cam,
            om = Omcc(:,pp);
            T = Tcc(:,pp);
            hand = handcc(pp);
            center2d(:,:,pp) = project_points_mirror2(center3d,om,T,hand,fc,cc,zeros(5,1),alpha_c);
        end;
        center = compute_structure2(center2d,Omcc,Tcc,handcc,fc,cc,zeros(5,1),alpha_c,MaxIter);
    else
        for pp = 1:n_cam,
            fc = fc_mat(:,pp);
            cc = cc_mat(:,pp);
            alpha_c = alpha_vec(pp);
            om = Omcc(:,pp);
            T = Tcc(:,pp);
            hand = handcc(pp);
            center2d(:,:,pp) = project_points_mirror2(center3d,om,T,hand,fc,cc,zeros(5,1),alpha_c);
        end;
        center = compute_structure2(center2d,Omcc,Tcc,handcc,fc_mat,cc_mat,zeros(5,n_cam),alpha_vec,MaxIter);
    end;
    err_X = center-center3d;
    std_X = std(err_X(:,~isnan(err_X(1,:))),0,2);
end;

save_name = [imgdir '/kinematic_animation.mat'];
if exist(save_name,'file')==2,
    fprintf(1,'\nLoading kinematics parameters file ''kinematic_animation.mat'' ...\n');
    load(save_name);
else
    fprintf(1,'\nERROR: ''kinematic_animation.mat'' not found!\n');
    return;
end;

filepath = [imgdir '/old'];
if ~exist(filepath,'dir'),
    return;
end;

ind_activeN = ind_active;
ind_downN = ind_down;
ind_upN = ind_up;
ind_maxdownN = ind_maxdown;
ind_maxupN = ind_maxup;
axes_meanN = axes_mean;
member_centerN = member_center;
member_axesN = member_axes;
member_omegaN = member_omega;
member_eulerAnglesN = member_eulerAngles;
ctpositionN = ctposition;
strphiN = strphi;
strthetaN = strtheta;
strpsiN = strpsi;
yawN = yaw;
pitchN = pitch;
rollN = roll;
wing_rootN = wing_root;
RbsN = Rbs;
R0N = R0;
body2strpN = body2strp;
strp_angleN = strp_angle;

save_name = [filepath '/kinematic_animation.mat'];
if exist(save_name,'file')==2,
    fprintf(1,'\nLoading original kinematics parameters file ...\n');
    load(save_name);
else
    fprintf(1,'\nERROR: Original ''kinematic_animation.mat'' not found in directory ''old''!\n');
    return;
end;

assert(isequal(ind_active,ind_activeN), 'The active frames of reprojection must be same with original!');
if ~isequal(ind_down,ind_downN),
    warning('Timing of downstroke is not same as original!');
end;
if ~isequal(ind_up,ind_upN),
    warning('Timing of upstroke is not same as original!');
end;
if ~isequal(ind_maxdown,ind_maxdownN),
    warning('Timing of maximum span in downstroke is not same as original!');
end;
if ~isequal(ind_maxup,ind_maxupN),
    warning('Timing of maximum span in upstroke is not same as original!');
end;

err_axes_mean = axes_meanN-axes_mean;
err_center = member_centerN(:,:,ind_active)-member_center(:,:,ind_active);
std_center = std(err_center,0,3);
err_axes = member_axesN(:,:,ind_active)-member_axes(:,:,ind_active);
err_omega = member_omegaN(:,:,ind_active)-member_omega(:,:,ind_active);
err_eulerAngles = member_eulerAnglesN(:,:,ind_active)-member_eulerAngles(:,:,ind_active);
std_eulerAngles = std(err_eulerAngles,0,3);
Q0 = trans_quat_axis(reshape(member_omega(:,:,ind_active), 3,[]));
Qt = trans_quat_axis(reshape(member_omegaN(:,:,ind_active), 3,[]));
% ind = sum(Q0.*Qt,1)<0;
% Qt(:,ind) = -Qt(:,ind);
% Q0(1:3,) = -Q0(1:3,:) 
% Qt = quatmul(Qt, Q0, 1);
% err_orientation = 2*acosd(Qt(4,:));
err_orientation = reshape(2*acosd(abs(sum(Q0.*Qt,1))),1,nparts+1,[]);
std_orientation = std(err_orientation,0,3);
err_wroot = wing_rootN-wing_root;
Qt = trans_quat_mat(RbsN*Rbs');
err_Rbs = 2*acosd(Qt(4));
Qt = trans_quat_mat(R0N*R0');
err_R0 = 2*acosd(Qt(4));
err_b2s = body2strpN-body2strp;
err_strp = strp_angleN-strp_angle;
std_strp = std(err_strp,0,2);

% show kinematic errors of body
nx = 1280;
ny = 800;
height = 4;

n = 6;
a = [ts(1),ts(end)]*1e3;
d = 0.1;    % margin of figure
w = 1-d*1.5;    % width and height of figure
h = (w-(n-1)*d/2)/n;
c = h+d/2;
c = (0:n-1)*c+d;
w1 = 0.5;       % line width of original data
w2 = 0.5;       % line width of reprojection data
e = a(1)-(a(2)-a(1))*d/(2*w);   % x postion of ylabel
palette = lines(7);
m = 12;     % font size

figure;
XX = [R0'*reshape(ctposition(:,1,:),3,[]); R0N'*reshape(ctpositionN(:,1,:),3,[])];
XX = XX-XX(:,1)*ones(1,nn);
b = [min(XX(:)), max(XX(:))];
delta = (b(2)-b(1))/10;
b = round(b+[-1,1]*delta);
y = 'ZYX';
for ii=1:3,
    axes('position',[d c(ii) w h]);
    jj = 3-ii+1;
    plot(1e3*ts, XX(jj,:), '-', 'color', palette(jj,:), 'linewidth',w1);   % original
    hold on;
    plot(1e3*ts, XX(jj+3,:), '-', 'color', palette(7,:), 'linewidth',w2);   % reprojected
    ylabel(['\fontname{Arial}\fontsize{' num2str(m) '}' y(ii) '_b (mm)'], 'position',[e, sum(b)/2]);
    if ii==1,
        xlabel(['\fontname{Arial}\fontsize{' num2str(m) '}time (ms)']);
        set(gca,'YMinorTick','on','box','off','xlim',a,'ylim',b);
    else
        set(gca,'XTick',[],'xcolor','w','YMinorTick','on','box','off','xlim',a,'ylim',b);
    end;
end;

XX = [yaw; pitch; roll; yawN; pitchN; rollN];
b = [min(XX(:)), max(XX(:))];
delta = (b(2)-b(1))/10;
b = round(b+[-1,1]*delta);
y = {'roll','pitch','yaw'};
for ii=1:3,
    axes('position',[d c(ii+3) w h]);
    jj = 3-ii+1;
    plot(1e3*ts, XX(jj,:), '-', 'color', palette(jj,:), 'linewidth',w1);   % original
    hold on;
    plot(1e3*ts, XX(jj+3,:), '-', 'color', palette(7,:), 'linewidth',w2);   % reprojected
    ylabel(['\fontname{Arial}\fontsize{' num2str(m) '}' y{ii} ' (^\circ)'], 'position',[e, sum(b)/2]);
    set(gca,'XTick',[],'xcolor','w','YMinorTick','on','box','off','xlim',a,'ylim',b);
end;
resolution = round(ny/height);   % dpi
set(gcf,'color','w','PaperPositionMode','Auto', 'PaperUnits','inches','PaperPosition',[0 0 nx ny]*2/resolution);
drawnow;
save_name = [imgdir '/kinematics_errb'];
print(gcf,'-djpeg',['-r' num2str(resolution*2)],save_name);
saveas(gcf, save_name,'fig');

% show kinematic errors of wings
figure;
XX = [strpsi; strpsiN];
b = [min(XX(:)), max(XX(:))];
delta = (b(2)-b(1))/10;
b = round(b+[-1,1]*delta);
y = 'lr';
for ii=1:2,
    axes('position',[d c(ii) w h]);
    jj = 2-ii+1;
    plot(1e3*ts, XX(jj,:), '-', 'color', palette(jj+1,:), 'linewidth',w1);   % original
    hold on;
    plot(1e3*ts, XX(jj+2,:), '-', 'color', palette(7,:), 'linewidth',w2);   % reprojected
    ylabel(['\fontname{Arial}\fontsize{' num2str(m) '}\psi_' y(ii) ' (^\circ)'],'position',[e, sum(b)/2]);
     if ii==1,
        xlabel(['\fontname{Arial}\fontsize{' num2str(m) '}time (ms)']);
        set(gca,'YMinorTick','on','box','off','xlim',a,'ylim',b);
    else
        set(gca,'XTick',[],'xcolor','w','YMinorTick','on','box','off','xlim',a,'ylim',b);
    end;
end;

XX = [strtheta; strthetaN];
b = [min(XX(:)), max(XX(:))];
delta = (b(2)-b(1))/10;
b = round(b+[-1,1]*delta);
for ii=1:2,
    axes('position',[d c(ii+2) w h]);
    jj = 2-ii+1;
    plot(1e3*ts, XX(jj,:), '-', 'color', palette(jj+1,:), 'linewidth',w1);   % original
    hold on;
    plot(1e3*ts, XX(jj+2,:), '-', 'color', palette(7,:), 'linewidth',w2);   % reprojected
    ylabel(['\fontname{Arial}\fontsize{' num2str(m) '}\theta_' y(ii) ' (^\circ)'], 'position',[e, sum(b)/2]);
    set(gca,'XTick',[],'xcolor','w','YMinorTick','on','box','off','xlim',a,'ylim',b);
end;

XX = [strphi; strphiN];
b = [min(XX(:)), max(XX(:))];
delta = (b(2)-b(1))/10;
b = round(b+[-1,1]*delta);
for ii=1:2,
    axes('position',[d c(ii+4) w h]);
    jj = 2-ii+1;
    plot(1e3*ts, XX(jj,:), '-', 'color', palette(jj+1,:), 'linewidth',w1);   % original
    hold on;
    plot(1e3*ts, XX(jj+2,:), '-', 'color', palette(7,:), 'linewidth',w2);   % reprojected
    ylabel(['\fontname{Arial}\fontsize{' num2str(m) '}\phi_' y(ii) ' (^\circ)'], 'position',[e, sum(b)/2]);
    set(gca,'XTick',[],'xcolor','w','YMinorTick','on','box','off','xlim',a,'ylim',b);
end;

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
save_name = [imgdir '/kinematics_errw'];
print(gcf,'-djpeg',['-r' num2str(resolution*2)],save_name);
saveas(gcf, save_name,'fig');

fprintf(1,'\nSaving difference of reprojection and original in ''kinematic_error.mat''.\n');
save_name = [imgdir '/kinematic_error.mat'];
string_save = ['save ' save_name ' n_frame N_active ind_active FramePS cyclePS FramePC Nstroke ' ...
                    'nn ts dt ind_down ind_up ind_maxdown ind_maxup axes_mean member_center ' ...
                    'member_axes member_omega ctposition member_eulerAngles wing_root Rbs R0 ' ...
                    'body2strp strp_angle strphi strtheta strpsi yaw pitch roll ind_downN ind_upN ' ...
                    'ind_maxdownN ind_maxupN axes_meanN member_centerN member_axesN member_omegaN ' ...
                    'member_eulerAnglesN ctpositionN strphiN strthetaN strpsiN yawN pitchN rollN ' ...
                    'wing_rootN RbsN R0N body2strpN strp_angleN err_axes_mean err_center std_center ' ...
                    'err_axes err_omega err_eulerAngles std_eulerAngles err_orientation std_orientation ' ...
                    'err_wroot err_Rbs err_R0 err_b2s err_strp std_strp'];
eval(string_save);
fprintf(1,'done...\n');
