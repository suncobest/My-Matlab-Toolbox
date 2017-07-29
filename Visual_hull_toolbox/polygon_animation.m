% find the codedir of the visual hull toolbox
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

save_name = [imgdir '/kinematic_animation.mat'];
if exist(save_name,'file')==2,
    fprintf(1,'\nLoading kinematics parameters file ''kinematic_animation.mat'' ...\n');
    load(save_name);
else
    fprintf(1,'\nERROR: ''kinematic_animation.mat'' not found!\n');
    return;
end;

% check and load stl model
load_stl = 1;
if ~exist('XYZb', 'var') || ~exist('XYZw','var') || isempty(XYZb) || isempty(XYZw),
    load_stl = 0;
else
    fprintf(1,'\nPolygon variable of body and wing detected!\n');
    flag = input('Do you want to change the model files or not? ([]=no, other=yes) ','s');
    if ~isempty(flag),
        fprintf(1,'\nYou choose to load another STL model ...\n');
        load_stl = 0;
    end;
end;

if ~load_stl,
    fprintf(1,'\nPlease make sure the STL models of body and wing are in the directory of ''stlTools''!\n');
    flag = input('Have you prepared STL models yet? ([]=yes, other=no) ','s');
    if ~isempty(flag),
        fprintf(1,'\nSTL models of body and wing are needed afterward...\n');
        return;
    end;
    % find the codedir of the stlTools toolbox
    codeDir = which('stlRead');
    if isempty(codeDir),
        fprintf(1,'\nYou must add the toolbox ''stlTools'' in your MATLAB path!\n');
        return;
    else
        flag = find(codeDir=='\' | codeDir=='/',1,'last');
        if ~isempty(flag),
            codeDir = codeDir(1 : flag(end)-1);
        end;
        fprintf(1,'\n\nDisplay all STL files in path of stlTools:\n');
        dir([codeDir '/*.stl']);
        if isempty(dir([codeDir '/*.stl'])),
            fprintf(1,'\nNo STL model was found in the path of ''stlTools''!\n');
            return;
        end;
    end;
    
    fprintf(1,'\nPlease input which type of model nomalization to perform ...\n');
    flag = input('Nomalize along the long axis or all 3 axes? ([]=long axis, other=all 3 axes) ','s');
    norm_sw = isempty(flag);
end;

while ~load_stl,
    save_name= input('Name of body model (without stl suffix): ( []=''bee_body'' ) ','s');
    if isempty(save_name),
        save_name = 'bee_body';
    end;
    if ~exist([codeDir '/' save_name '.stl'],'file'),
        fprintf(1,'\nNo file %s found in the path! Please input again!\n', save_name);
        continue;
    end;
    % load and normalize model of body and wing(ascii or binary STL)
    [XYZb,faceb] = stlRead([codeDir '/' save_name '.stl']);
    if isempty(XYZb) || isempty(faceb),
        fprintf(1,'\nSTL file of body found empty! Please input again!\n');
        continue;
    end;
    nptb = size(XYZb,1);
    if norm_sw,
        XYZb =2*(XYZb/(max(XYZb(:,1))-min(XYZb(:,1))))';         % normalized by body length (radius to 1)
    else
        XYZb = 2*(XYZb./(ones(nptb,1)*(max(XYZb,[],1)-min(XYZb, [],1))))';         % normalized like cube (radius to 1)
    end;
    
    save_name= input('Name of wing model (without stl suffix): ( []=''bee_wingr'' ) ','s');
    if isempty(save_name),
        save_name = 'bee_wingr';
    end;
    if ~exist([codeDir '/' save_name '.stl'],'file'),
        fprintf(1,'\nNo file %s found in the path! Please input again!\n', save_name);
        continue;
    end;
    [XYZw,facew] = stlRead([codeDir '/' save_name '.stl']);
    if isempty(XYZw) || isempty(facew),
        fprintf(1,'\nSTL file of wing found empty! Please input again!\n');
        continue;
    end;
    nptw = size(XYZw,1);
    if norm_sw,
        %     XYZw =2*(XYZw/max(XYZw(:,2)))';        % normalized by wing length from origin to tip (radius to 1)
        XYZw =2*(XYZw/(max(XYZw(:,2))-min(XYZw(:,2))))';        % normalized by wing length (radius to 1)
    else
        XYZw = XYZw(:,1:2);     % thickness of wing is set to zero
        XYZw = 2*[(XYZw./(ones(nptw,1)*(max(XYZw,[],1)-min(XYZw, [],1))))'; zeros(1,nptw)];        % normalized like square (radius to 1)
    end;
    load_stl = 1;
end;

% some settings
% az = -30;
% el = 30;
palette = lines(7);

nlen_vec = input('Length factor of principle axis vector to show: ([]=1.5 times semi_axes) ');
if isempty(nlen_vec),
    nlen_vec = 1.5;
end;

% make wings' roots symmetrical
wroot = wing_root;
wroot(2,2) = -wroot(2,2);
wroot = sum(wroot,2)/2*ones(1,2);
wroot(2,2) = -wroot(2,2);

% model animation
h = figure;
for kk=1:nn,
    figure(h);hold off;
    % plot body
    ax = axes_mean(:,1);
    ct0 = ctposition(:,1,kk);
    vc0 = rodrigues(axisAngle(:,1,kk));
    % plot axes
    ctt = ct0(:,ones(1,2));
    xyz = ctt+vc0(:,1:2)*diag(ax(1:2))*nlen_vec;
    for i=1:2, 
        plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',2);
        hold on;
    end;
    % resize, reorient, and relocate polygon
    if norm_sw,
        XYZe = vc0*ax(1)*XYZb+ct0(:,ones(1,nptb));
    else
        XYZe = vc0*diag(ax)*XYZb+ct0(:,ones(1,nptb));
    end;
    patch('Faces', faceb, 'Vertices', [XYZe(1,:); XYZe(3,:); -XYZe(2,:)]','FaceColor',palette(1,:),...
        'EdgeColor',palette(1,:),'FaceAlpha',0.1);
    
    % plot right wing
    ax = axes_mean(:,2);
    ct = ct0+vc0*wroot(:,1);
    vc = rodrigues(axisAngle(:,2,kk));
    % plot axes
    ctt = ct(:,ones(1,2));
    xyz = ctt+vc(:,1:2)*diag(ax(1:2))*nlen_vec;
    for i=1:2,
        plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',2);
    end;
    % resize, reorient, and relocate polygon
    if norm_sw,
        XYZe = vc*ax(2)*XYZw+ct(:,ones(1,nptw));
    else
        XYZe = vc*diag(ax)*XYZw+ct(:,ones(1,nptw));
    end;
    patch('Faces', facew, 'Vertices', [XYZe(1,:); XYZe(3,:); -XYZe(2,:)]','FaceColor',palette(2,:),...
        'EdgeColor',palette(2,:),'FaceAlpha',0.1);
    
    % plot left wing
    ax = axes_mean(:,3);
    ax(2) = -ax(2);     % mirror y axis
    ct = ct0+vc0*wroot(:,2);
    vc = rodrigues(axisAngle(:,3,kk));
    % plot axes
    ctt = ct(:,ones(1,2));
    xyz = ctt+vc(:,1:2)*diag(ax(1:2))*nlen_vec;
    for i=1:2,
        plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',2);
    end;
    
    % resize, reorient, and relocate polygon
    if norm_sw,
        XYZe = vc*ax(2)*diag([-1,1,-1])*XYZw+ct(:,ones(1,nptw));
    else
        XYZe = vc*diag(ax)*XYZw+ct(:,ones(1,nptw));
    end;
    patch('Faces', facew, 'Vertices', [XYZe(1,:); XYZe(3,:); -XYZe(2,:)]','FaceColor',palette(3,:),...
        'EdgeColor',palette(3,:),'FaceAlpha',0.1);
    
    axis equal; grid on;
    set(gcf, 'renderer', 'zbuffer','color',[1,1,1]*0.7);
    cameratoolbar('ResetCameraAndSceneLight');
    cameratoolbar('Togglescenelight');
    view(az,el);
    axis(axes_range);
    pause(0.01);
end;
disp('Done...');
