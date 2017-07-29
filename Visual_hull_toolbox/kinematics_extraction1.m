% This script will redefine orientation of wings (swap x and y axes)
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
    fprintf(1,'\nERROR: ''segment_variables.mat'' not found!\n');
    return;
end;

% view axis for 3D space
nlen_vec = input('Length factor of principle axis vector to show: ([]=1.5 times semi_axes) ');
if isempty(nlen_vec),
    nlen_vec = 1.5;
end;

ind = ~any(isnan(member_nax));
Nb = sum(ind);
assert(isequal(find(ind),ind_active) && Nb==n_frame, 'All image frames must be active to extract kinematics!');
semiax = member_axes./(repmat(reshape(member_nax,1,nparts+1,[]),[3,1]));
axes_mean = sum(semiax,3)/Nb;
axis_ux = NaN(3,nparts+1,n_frame);
axis_uy = axis_ux;
for kk=1:n_frame,
    for ii=1:nparts+1,
        R = rodrigues(member_omega(:,ii,kk));
        axis_ux(:,ii,kk) = R(:,1);
        axis_uy(:,ii,kk) = R(:,2);
    end;
end;
meanhead = member_center + repmat(axes_mean(1,:),[3,1,n_frame]).*axis_ux;
meanfoot = member_center - repmat(axes_mean(1,:),[3,1,n_frame]).*axis_ux;

flag = 1;
while flag,
    FramePS = input('What is the frame rate (FPS) of the image sequence? ');
    if isempty(FramePS) || FramePS<=0,
        fprintf(1,'\nUnexpected input!\n');
        continue;
    end;
    flag = 0 ;
end;


if nparts==2,
    flag = 1;
    while flag,
        cyclePS = input('What is the motion frequency (cycle per second) of the animal? ');
        if isempty(cyclePS) || cyclePS<=0,
            fprintf(1,'\nUnexpected input!\n');
            continue;
        end;
        flag = 0 ;
    end;
    FramePC = FramePS/cyclePS;
    assert(FramePC>=8,'Too few images in one stroke, you shoud raise the shooting speed!');
    flag = 1;
    while flag,
        idk = input('Frame number at the start of the 1st downstroke: ( []=1 ) ');
        if isempty(idk),
            idk = 1;
        else
            idk = round(idk);
            if idk<1 || idk>n_frame,
                fprintf(1,'\nUnexpected input!\n');
                continue;
            end;
        end;
        flag = 0;
    end;
    ind_reversal = round(idk : FramePC/2 : n_frame);
    assert(length(ind_reversal)>=3,'There must be at least one whole stroke!');
    % load the last frame of point cloud
    count = ceil(n_frame/nfpu);
    base = (count-1)*nfpu;
    nc = sprintf(['%0' ndigit 'd'],count);
    save_name = [imgdir '/segment_part' nc '.mat'];
    if exist(save_name,'file')==2,
        load(save_name);
    else
        fprintf(1,'\nERROR: cannot find data ''%s''!\n',save_name);
        return;
    end;
    bk = n_frame-base;
    XX = X_cell{bk};
    seg_id = seg_cell{bk};
    
    figure(3); hold off;
    idx = (seg_id==0);
    ii = nparts+2;
    plot3(XX(1,idx),XX(3,idx),-XX(2,idx), '.', 'color',palette(ii,:));
    hold on;
    for ii=1:nparts+1,
        idx = (seg_id==ii);
        plot3(XX(1,idx),XX(3,idx),-XX(2,idx), '.', 'color',palette(ii,:));
    end;
    axis equal vis3d off;
    view(az,el);

    % refine body axes: compute lateral direction of body (from left to right wing joint)
    body_y = reshape(meanfoot(:,2,:)-meanfoot(:,3,:), 3,[]);
    body_y = body_y./repmat(sqrt(sum(body_y.*body_y)),3,1);    % axis2 of body
    body_x = reshape(axis_ux(:,1,:), 3,[]);   % axis1 of body
    fprintf(1,['\nRefine the body lateral orientation if there are two members (wings) besides body:\n'...
        'body is in the middle (blue), seperating two wings (green and red).\n']);
    flag = input('Please identify the color of right wing: ([]=green, other=red) ','s');
    if ~isempty(flag),      % right wing: 2nd column, left wing: 3rd column
        body_y = -body_y;
        member_nax = member_nax([1 3 2],:);
        member_center = member_center(:,[1 3 2],:);
        member_axes = member_axes(:,[1 3 2],:);
        member_omega = member_omega(:,[1 3 2],:);
        meanhead = meanhead(:,[1 3 2],:);
        meanfoot = meanfoot(:,[1 3 2],:);
        axes_mean = axes_mean(:,[1 3 2]);
        axis_ux = axis_ux(:,[1 3 2],:);
        axis_uy = axis_uy(:,[1 3 2],:);
    end;
    % using middle downstroke to denoise body's lateral direction
    flag = 1;
    while flag,
        idk = input('Frame number with greatest span during the 1st downstroke: ');
        idk = round(idk);
        if  length(idk)~=1|| idk<1 || idk>n_frame,
            fprintf(1,'\nUnexpected input!\n');
            continue;
        end;
        flag = 0;
    end;
    ind_maxspan = round(idk(1) : FramePC : n_frame);
    Nb = length(ind_maxspan);
    if FramePC>10,
        ind = [ind_maxspan-1; ind_maxspan; ind_maxspan+1];
        a = zeros(3,Nb,3);
        b = a;
        a(:,:,2) = body_x(:,ind_maxspan);
        b(:,:,2) = body_y(:,ind_maxspan);
        idx = ind(1,:);
        if idx(1)<1,
            idx(1) = 1;
        end;
        a(:,:,1) = body_x(:,idx);
        b(:,:,1) = body_y(:,idx);
        idx = ind(3,:);
        if idx(end)>n_frame,
            idx(end)=n_frame;
        end;
        a(:,:,3) = body_x(:,idx);
        b(:,:,3) = body_y(:,idx);
        % mean of three orientation
        a = mean(a,3);
        b = mean(b,3);
        a = a./repmat(sqrt(sum(a.*a)),3,1);
        % orthogonalize the lateral direction with body's long axis
        b = b-repmat(sum(b.*a),3,1).*a;
        b = b./repmat(sqrt(sum(b.*b)),3,1);
    else
        %  frame numbers of middle downstroke and upstroke
        a = body_x(:,ind_maxspan);
        b = body_y(:,ind_maxspan);
        b = b-repmat(sum(b.*a),3,1).*a;
        b = b./repmat(sqrt(sum(b.*b)),3,1);
    end;
    ind = ind_maxspan;
    omega = zeros(3,Nb);
    for ii = 1:Nb,
        omega(:,ii) = rodrigues([a(:,ii), b(:,ii), cross(a(:,ii),b(:,ii))]);
    end;
    % quaternion interpolation (squad function) with clamped orientation
    if ind(1)>1,
        ind = [1,ind];
        omega = [omega(:,1),omega];
    end;
    if ind(end)<n_frame,
        ind = [ind, n_frame];
        omega = [omega, omega(:,end)];
    end;
    omega = trans_quat_axis(squad(ind, trans_quat_axis(omega),1:n_frame));
    % Re-orthogonalize the lateral direction with original body main axis
    flag = input('Call back original body''s main axis or not? ([]=no (recommended), other=yes) ','s');
    if ~isempty(flag),
        body_y = zeros(3,n_frame);
        for kk=1:n_frame,
            ax = rodrigues(omega(:,kk));
            body_y(:,kk) = ax(:,2);
        end;
        body_y = body_y-repmat(sum(body_y.*body_x),3,1).*body_x;
        body_y = body_y./repmat(sqrt(sum(body_y.*body_y)),3,1);
        for kk = 1:n_frame,
            a = body_x(:,kk);
            b = body_y(:,kk);
            omega(:,kk) = rodrigues([a, b, cross(a,b)]);
        end;
    end;
    % update the body orientation
    member_omega(:,1,:) = omega;   
    
    % compute wings' center, tip, root position and orientation with respect to body coordinate system
    temp = repmat(member_center(:,1,:),[1, 2]);
    wcenter = member_center(:,2:3,:)-temp;
    whead = meanhead(:,2:3,:)-temp;
    wfoot = meanfoot(:,2:3,:)-temp;
    temp = axes_mean(1,1);
    XX = [wcenter,whead,wfoot,axis_ux(:,2:3,:),axis_uy(:,2:3,:)];
    for kk = 1:n_frame,
        ax = rodrigues(omega(:,kk));
        a = ax(:,1);
        b = ax(:,2);
        body_x(:,kk) = a;
        body_y(:,kk) = b;
        axis_ux(:,1,kk) = a;
        axis_uy(:,1,kk) = b;
        meanhead(:,1,kk) = member_center(:,1,kk) + temp*a;
        meanfoot(:,1,kk) = member_center(:,1,kk) - temp*a;
        XX(:,:,kk) = ax'*XX(:,:,kk);
    end;
    wcenter = reshape(XX(:,1:2,:),3,[]);
    whead = XX(:,3:4,:);
    wfoot = XX(:,5:6,:);
    wing_x = reshape(XX(:,7:8,:),3,[]);
    wing_y = reshape(XX(:,9:10,:),3,[]);

    % compute wing joints
    [rw_root, rlen_rc, ~, rw_std] = lines_joint(wcenter(:,1:2:end), wing_x(:,1:2:end));
    [lw_root, llen_rc, ~, lw_std] = lines_joint(wcenter(:,2:2:end), wing_x(:,2:2:end));
    rlen_rc = mean(rlen_rc,2);
    llen_rc = mean(llen_rc,2);
    wing_root = mean(wfoot,3);
    
    % draw wing root and tip in body's coordinate system
    a = reshape(whead,3,[]);
    b = reshape(wfoot,3,[]);
    wing_root1 = [rw_root,lw_root];
    figure(2); hold off;
    arrow3(temp*[0 0 -1],temp*[0 0 1],'b2',20,4);
    hold on;
    plot3([a(2,:);b(2,:)],[a(3,:);b(3,:)],[a(1,:);b(1,:)]);
    plot3(wing_root(2,:),wing_root(3,:),wing_root(1,:),'go','markersize',10);
    plot3(wing_root1(2,:),wing_root1(3,:),wing_root1(1,:),'r+','markersize',10);
    axis image; grid on;
    view(az,el);
    camlight('headlight');
    lighting gouraud;
    
    % kinematics measurement
    % compute stroke plane (xy plane) in the body coordinate system (identity matrix)
    %     ind = ind_reversal(1:2);   % two extremes of the 1st stroke
    if mod(length(ind_reversal),2),   % two extremes of all strokes
        ind = ind_reversal(1:end-1);
    else
        ind = ind_reversal;
    end;
    a = reshape([whead(:,:,ind),wfoot(:,:,ind)],3,[]);
    a = [a; ones(1,size(a,2))];
    a = a*a';
    [~,~,b] = svd(a);
    b = b(1:3,4);
    a = b/norm(b);
    if a(1)>0,      % the vertical z axis of stroke plane must be at the minus half of body's x axis
        a = -a;
    end;
    b = [a(3); 0; -a(1)];  % cross([0;1;0],a);
    b = b/norm(b); % the x axis on stroke plane
    body2strp = -acosd(b(1));       % rotation around y axis of body (+: counter-clockwise)
    % wing's orientation wrt stroke plane coordinate system
    Rbs = rodrigues([0;1;0]*body2strp*pi/180);
    wing_x = Rbs'*wing_x;
    wing_y = Rbs'*wing_y;

    % compute orientation of wing wrt body (stroke plane): stroke angle, deviation angle and incidence
    b = [wing_x(1:2,:); zeros(1,n_frame*2)];
    b = b./(ones(3,1)*sqrt(sum(b.*b)));
    temp = b(2,:);
    temp(2:2:end) = -temp(2:2:end);
    stroke_phi = reshape(atan2d(b(1,:),temp),2,[]);   % stroke angle
    stroke_theta = -reshape(asind(wing_x(3,:)),2,[]);   % deviation angle
    % right wing: cross(wing_x(:,1:2:end), [0;0;1]*ones(1,n_frame));
    % left wing:  cross( [0;0;1]*ones(1,n_frame), wing_x(:,2:2:end));
    b = [wing_x(2,:); -wing_x(1,:); zeros(1,n_frame*2)];
    b(:,2:2:end) = -b(:,2:2:end);
    b = b./(ones(3,1)*sqrt(sum(b.*b)));
    c = cross(wing_x,b);
    c(:,2:2:end) = -c(:,2:2:end);
    stroke_psi = reshape(atan2d(sum(wing_y.*c),sum(wing_y.*b)),2,[]);   % incidence: acosd(sum(wing_y.*b));
    ind = stroke_psi<-90;
    stroke_psi(ind) = stroke_psi(ind)+360;   % incidence range: [-90, 270)
    
    % compute orientation of body
    c = [0;1;0]*ones(1,n_frame);    % the z axis of initial coordinate system
    a = [body_x(1,:); zeros(1,n_frame); body_x(3,:)];
    a = a./(ones(3,1)*sqrt(sum(a.*a)));     % the intermediate x axis
    b= cross(c(:,1),a(:,1));    % the y axis of initial coordinate system
    R0 = [a(:,1),b,c(:,1)];      % the initial coordinate system (laboratory-fixed)
    body_yaw = atan2d(sum(a.*(b*ones(1,n_frame))), sum(a.*(a(:,1)*ones(1,n_frame))));
    body_pitch = -asind(sum(body_x.*c));
    b = cross(c, body_x);
    b =  b./(ones(3,1)*sqrt(sum(b.*b)));     % the intermediate y axis
    c = cross(body_x, b);     % the intermediate z axis
    body_roll = atan2d(sum(body_y.*c), sum(body_y.*b));
    strp_angle = -body2strp-body_pitch;
    
    fracw = input('Fraction number to approximate thickness of wings: ([]=0.2) ');
    if isempty(fracw),
        fracw = 0.2;
    end;
    axes_mean(3,2:3) = axes_mean(3,2:3)*fracw;
    
    flag = 1;
    while flag,
        ndiv = input('Subdivision number for one stroke to interpolate animation: ([]=200) ');
        if isempty(ndiv),
            ndiv = 200;
        else
            ndiv = round(ndiv);
            if  ndiv<=40,
                fprintf(1,'\nSubdivision number for one stroke must not be less than 40!\n');
                continue;
            end;
        end;
        flag = 0;
    end;

    % spline interpolation
    Nstroke = (n_frame-1)/FramePC;
    t0 = (0:n_frame-1)/FramePS;        % unit: s
    nn = round(Nstroke*ndiv+1);
    ts = linspace(t0(1),t0(end),nn);
    dt = ts(2)-ts(1);
    XX = [stroke_phi; stroke_theta; stroke_psi; body_yaw; body_pitch; body_roll; ...
        reshape(member_center,3*(nparts+1),[])];
    [a, b, c] = spline_eval(t0, spline_interp3(t0, XX), ts);
    % orientation of wing wrt body (stroke plane)
    strphi1 = a(1:2,:);       % unit: degree
    strphi1dt = b(1:2,:);       % unit: degree/s
    strphi1dt2 = c(1:2,:);      % unit: degree/s^2
    strtheta1 = a(3:4,:);
    strtheta1dt = b(3:4,:);
    strtheta1dt2 = c(3:4,:);
    strpsi1 = a(5:6,:);
    strpsi1dt = b(5:6,:);
    strpsi1dt2 = c(5:6,:);
    % orientation of body
    yaw1 = a(7,:);
    yaw1dt = b(7,:);
    yaw1dt2 = c(7,:);
    pitch1 = a(8,:);
    pitch1dt = b(8,:);
    pitch1dt2 = c(8,:);
    roll1 = a(9,:);
    roll1dt = b(9,:);
    roll1dt2 = c(9,:);
    % position of member center
    ctposition = reshape(a(10:end,:),3,nparts+1,[]);
    ctspeed = reshape(b(10:end,:),3,nparts+1,[]);
    ctacceler = reshape(c(10:end,:),3,nparts+1,[]);
    
   %% recompute orientation kinematics after quaternion interplation
    % redefine wings' axes to suit Euler angles:
    % swap axis x and y for right wing; swap x and y and then reverse y axis for left wing
    axes_mean(1:2,2:3) = axes_mean([2 1],2:3);
    member_axes(1:2,2:3,:) = member_axes([2 1],2:3,:);
    for kk = 1:n_frame,
        R = rodrigues(member_omega(:,2,kk));        % adjust orientation of right wing
        member_omega(:,2,kk) = rodrigues(R*[0 1 0;1 0 0;0 0 -1]);
        R = rodrigues(member_omega(:,3,kk));        % adjust orientation of left wing
        member_omega(:,3,kk) = rodrigues(R*[0 -1 0;1 0 0;0 0 1]);
    end;
    temp = axis_uy(:,2:3,:);
    axis_uy(:,2:3,:) = axis_ux(:,2:3,:);
    axis_uy(:,3,:) = -axis_uy(:,3,:);
    axis_ux(:,2:3,:) = temp;
    temp = wing_y;
    wing_y = wing_x;
    wing_y(:,2:2:end) = -wing_y(:,2:2:end);
    wing_x = temp;
    % quaternion interpolation:
    % prefix 'ex' denotes extrinsic system, 'in' denotes intrinsic system
    axisAngle = NaN(3,nparts+1,nn);
    ex_rotspeed = axisAngle;
    in_rotspeed = axisAngle;
    XX = zeros(9,nparts+1,nn);
    for ii = 1:nparts+1,
        Q0 = trans_quat_axis(reshape(member_omega(:,ii,:),3,[]));
        Qt = squad(t0, Q0, ts);
        axisAngle(:,ii,:) =  trans_quat_axis(Qt);
        ax = zeros(9,nn);
        a = ax;
        for kk=1:nn,
            ax(:,kk) = reshape(trans_quat_mat(Qt(:,kk)),9,[]);
        end;
        XX(:,ii,:) = ax;
        % derivative of ax with time
        a(:,1) = (ax(:,2)-ax(:,1))/dt;
        a(:,end) = (ax(:,end)-ax(:,end-1))/dt;
        a(:,2:end-1) = (ax(:,3:end)-ax(:,1:end-2))/(2*dt);
        for kk = 1:nn,
            R = reshape(ax(:,kk),3,[]);
            temp = reshape(a(:,kk),3,[])*R';
            temp = (temp-temp')/2;
            temp = [temp(3,2); temp(1,3); temp(2,1)];
            ex_rotspeed(:,ii,kk) = temp;
            in_rotspeed(:,ii,kk) = R'*temp;
        end;
    end;
   
    body_xx = reshape(XX(1:3,1,:),3,[]);
    body_yy = reshape(XX(4:6,1,:),3,[]);
    wing_xx = reshape(XX(1:3,2:3,:),3,[]);
    wing_yy = reshape(XX(4:6,2:3,:),3,[]);
    for kk=1:nn,
        R = reshape(XX(:,1,kk),3,[]);
        ii = kk*2-1 : kk*2;
        wing_xx(:, ii) = R'*wing_xx(:, ii);
        wing_yy(:, ii) = R'*wing_yy(:, ii);
    end;
    wing_xx = Rbs'*wing_xx;
    wing_yy = Rbs'*wing_yy;
    
    % recompute orientation of wing wrt stroke plane
    b = [wing_yy(1:2,:); zeros(1,nn*2)];
    b = b./(ones(3,1)*sqrt(sum(b.*b)));
    temp = b(1,:);
    temp(2:2:end) = -temp(2:2:end);
    strphi = reshape(atan2d(temp,b(2,:)),2,[]);   % stroke angle
    temp = asind(wing_yy(3,:));
    temp(1:2:end) = -temp(1:2:end);
    strtheta = reshape(temp,2,[]);   % deviation angle
    b = [wing_yy(2,:); -wing_yy(1,:); zeros(1,nn*2)];
    b = b./(ones(3,1)*sqrt(sum(b.*b)));
    c = cross(wing_yy,b);
    strpsi = reshape(atan2d(sum(wing_xx.*c),sum(wing_xx.*b)),2,[]);   % incidence
    ind = strpsi<-90;
    strpsi(ind) = strpsi(ind)+360;   % incidence range: [-90, 270)
    
    % recompute orientation of body
    c = [0;1;0]*ones(1,nn);    % the z axis of initial coordinate system
    a = [body_xx(1,:); zeros(1,nn); body_xx(3,:)];
    a = a./(ones(3,1)*sqrt(sum(a.*a)));     % the intermediate x axis
    b= cross(c(:,1),a(:,1));    % the y axis of initial coordinate system
    yaw = atan2d(sum(a.*(b*ones(1,nn))), sum(a.*(a(:,1)*ones(1,nn))));
    pitch = -asind(sum(body_xx.*c));
    b = cross(c, body_xx);
    b =  b./(ones(3,1)*sqrt(sum(b.*b)));     % the intermediate y axis
    c = cross(body_xx, b);     % the intermediate z axis
    roll = atan2d(sum(body_yy.*c),sum(body_yy.*b));
    
    % derivative with time
    temp = [strphi; strtheta; strpsi; yaw; pitch; roll];
    XX = zeros(9,nn);
    XX(:,1) = (temp(:,2)-temp(:,1))/dt;
    XX(:,end) = (temp(:,end)-temp(:,end-1))/dt;
    XX(:,2:end-1) = (temp(:,3:end)-temp(:,1:end-2))/(2*dt);
    strphidt = XX(1:2, :);
    strthetadt = XX(3:4, :);
    strpsidt = XX(5:6, :);
    yawdt = XX(7, :);
    pitchdt = XX(8, :);
    rolldt = XX(9, :);
    
    % difference between spline interpolated Euler angle and quaternion interpolation
    ex_rotspeed1 = NaN(3,nparts+1,nn);
    XX = [yawdt; pitchdt; rolldt; yaw; pitch; roll; yaw1dt; pitch1dt; roll1dt; ...
        yaw1; pitch1; roll1; strphidt; strthetadt; strpsidt; strphi; strtheta; strpsi; ...
        strphi1dt; strtheta1dt; strpsi1dt; strphi1; strtheta1; strpsi1;]*pi/180;   %Unit: radian/s
    % change sign of the 1st and 2nd Euler angle of right wing due to the definition
    ind = [13, 15, 19, 21, 25, 27, 31, 33];
    XX(ind,:) = -XX(ind,:);
    % inav denotes intrinsic angular velocity, exav denotes extrinsic angular velocity;
    %Unit: radian/s. subscript b stands for body, l the left wing, r the right wing.
    [inavb, exavb] = deuler2rotspeed(XX(1:3,:),XX(4:6,:),'ZYX');
    [inavb1, exavb1] = deuler2rotspeed(XX(7:9,:),XX(10:12,:),'ZYX');
    inavr = deuler2rotspeed(XX(13:2:17,:),XX(19:2:23,:),'ZXY');
    inavl = deuler2rotspeed(XX(14:2:18,:),XX(20:2:24,:),'ZXY');
    inavr1 = deuler2rotspeed(XX(25:2:29,:),XX(31:2:35,:),'ZXY');
    inavl1 = deuler2rotspeed(XX(26:2:30,:),XX(32:2:36,:),'ZXY');
    exavb = R0*exavb;
    exavb1 = R0*exavb1;

    % error of Euler angle interpolation
    err_phi = std(strphi-strphi1,0,2);
    err_theta = std(strtheta-strtheta1,0,2);
    err_psi = std(strpsi-strpsi1,0,2);
    err_yaw = std(yaw-yaw1,0,2);
    err_pitch = std(pitch-pitch1,0,2);
    err_roll = std(roll-roll1,0,2);
    % error of change rate of Euler angle
    err_phidt = std(strphidt-strphi1dt,0,2);
    err_thetadt = std(strthetadt-strtheta1dt,0,2);
    err_psidt = std(strpsidt-strpsi1dt,0,2);
    err_yawdt = std(yawdt-yaw1dt,0,2);
    err_pitchdt = std(pitchdt-pitch1dt,0,2);
    err_rolldt = std(rolldt-roll1dt,0,2);
    % error of angular velocity wrt local axes in degree
    err_inavb = std(inavb-inavb1,0,2)*180/pi;
    err_inavr = std(inavr-inavr1,0,2)*180/pi;
    err_inavl = std(inavl-inavl1,0,2)*180/pi;
    
    % show orientation kinematics
    show_kinematics;
    
    % timing of 1st downstoke (right wing):
    [~,idk] = min(strphi(1,1:ndiv));
    Nb = floor((nn-idk)/ndiv);
    XX = [strphi; strtheta; strpsi];
    XX = reshape(XX(:, idk : idk+Nb*ndiv-1), 6, ndiv, []);
    a = mean(XX,3);
    b = sqrt(sum((XX-a(:,:,ones(1,Nb))).^2, 3)/(Nb-1));
    avphi = a(1:2,:);
    avtheta = a(3:4,:);
    avpsi = a(5:6,:);
    stdphi = b(1:2,:);
    stdtheta = b(3:4,:);
    stdpsi = b(5:6,:);
    
    [a, b] = max(XX(1:2,:,:),[],2);
    maxphi = mean(a,3);
    Tmaxphi = (mean(b,3)-1)/ndiv;
    minphi = mean(XX(1:2,1,:),3);
    % show wings' kinematics displayed in one period
    show_euler_average;
    
    
else
    
    % compute linear and angular velocity of every member for animation
    flag = 1;
    while flag,
        ndiv = input('Subdivision number for one frame to interpolate animation: ([]=10) ');
        if isempty(ndiv),
            ndiv = 10;
        else
            ndiv = round(ndiv);
            if  ndiv<=1,
                fprintf(1,'\nUnexpected input!\n');
                continue;
            end;
        end;
        flag = 0;
    end;
    t0 = (0:n_frame-1)/FramePS;        % unit: s
    nn = (n_frame-1)*ndiv+1;
    ts = linspace(t0(1),t0(end),nn);
    dt = ts(2)-ts(1);
    axisAngle = NaN(3,nparts+1,nn);
    ex_rotspeed = axisAngle;
    in_rotspeed = axisAngle;
    
    temp = spline_interp3(t0, reshape(member_center,3*(nparts+1),[]));
    [ctposition,ctspeed,ctacceler] = spline_eval(t0, temp, ts);
    ctposition = reshape(ctposition,3,nparts+1,[]);
    ctspeed = reshape(ctspeed,3,nparts+1,[]);
    ctacceler = reshape(ctacceler,3,nparts+1,[]);
    for ii = 1:nparts+1,
        Q0 = trans_quat_axis(reshape(member_omega(:,ii,:),3,[]));
        Qt = squad(t0, Q0, ts);
        axisAngle(:,ii,:) =  trans_quat_axis(Qt);
        ax = zeros(9,nn);
        a = ax;
        for kk=1:nn,
            ax(:,kk) = reshape(trans_quat_mat(Qt(:,kk)),9,[]);
        end;
        % derivative of ax with time
        a(:,1) = (ax(:,2)-ax(:,1))/dt;
        a(:,end) = (ax(:,end)-ax(:,end-1))/dt;
        a(:,2:end-1) = (ax(:,3:end)-ax(:,1:end-2))/(2*dt);
        for kk = 1:nn,
            R = reshape(ax(:,kk),3,[]);
            temp = reshape(a(:,kk),3,[])*R';
            temp = (temp-temp')/2;
            temp = [temp(3,2); temp(1,3); temp(2,1)];
            ex_rotspeed(:,ii,kk) = temp;
            in_rotspeed(:,ii,kk) = R'*temp;
        end;
    end;
end;

% mesh of unit sphere
Nb = 10;
npts = (Nb+1)^2;
[Xe, Ye, Ze] = sphere(Nb);
[faces,XYZs] = surf2patch(Xe,Ye,Ze);
XYZs = XYZs';
% show animation
a = min(reshape(member_center(:,1,:),3,[]),[],2)-axes_mean(1,1)*2;
b = max(reshape(member_center(:,1,:),3,[]),[],2)+axes_mean(1,1)*2;
axes_range = [a(1),b(1),a(3),b(3),-b(2),-a(2)];
ellipsoid_animation;


fprintf(1,'\nSaving kinematics parameters in ''kinematic_animation.mat''.\n');
save_name = [imgdir '/kinematic_animation.mat'];
string_save = ['save ' save_name ' n_frame ind_active active_images ndigit nx ny height palette color_cell axes_mean ' ...
    'nparts member_nax member_center member_axes member_omega axis_ux axis_uy meanhead meanfoot nfpu n_unit '...
    'FramePS ndiv nn ts dt ctposition ctspeed ctacceler axisAngle ex_rotspeed in_rotspeed axes_range az el'];
if nparts==2,
    string_save = [string_save ' cyclePS FramePC Nstroke body_x body_y ind_reversal ind_maxspan wcenter ' ...
        'whead wfoot wing_x wing_y wing_root wing_root1 Rbs R0 body2strp stroke_phi stroke_theta stroke_psi ' ...
        'body_yaw body_pitch body_roll strp_angle body_xx body_yy wing_xx wing_yy strphi strtheta strpsi ' ...
        'yaw pitch roll strphidt strthetadt strpsidt yawdt pitchdt rolldt strphi1 strtheta1 strpsi1 yaw1 pitch1 roll1 '...
        'strphi1dt strtheta1dt strpsi1dt yaw1dt pitch1dt roll1dt strphi1dt2 strtheta1dt2 strpsi1dt2 yaw1dt2 '...
        'pitch1dt2 roll1dt2 inavb inavr inavl inavb1 inavr1 inavl1 exavb exavb1 err_phi err_theta err_psi err_yaw '...
        'err_pitch err_roll err_phidt err_thetadt err_psidt err_yawdt err_pitchdt err_rolldt err_inavb err_inavr err_inavl '...
        'avphi avtheta avpsi stdphi stdtheta stdpsi maxphi minphi Tmaxphi'];
end;
eval(string_save);
fprintf(1,'done...\n');

