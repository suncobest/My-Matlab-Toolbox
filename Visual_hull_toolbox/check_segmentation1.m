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

nlen_vec = input('Length factor of principle axis vector to show: ([]=1.5 times semi_axes) ');
if isempty(nlen_vec),
    nlen_vec = 1.5;
end;

check_direction = input('Check direction of every member? ([]=yes, other=no) ','s');
check_direction = isempty(check_direction);

Nstart = input('Please set the start frame: ([]=1) ');
if isempty(Nstart),
    Nstart = 1;
end;
kk = Nstart;
count = ceil(kk/nfpu);
base = (count-1)*nfpu;
nc = sprintf(['%0' ndigit 'd'],count);
save_name = [imgdir '/segment_part' nc '.mat'];
if exist(save_name,'file')==2,
    load(save_name);
else
    fprintf(1,'\nERROR: cannot find data ''%s''!\n',save_name);
    return;
end;

for kk = Nstart:n_frame,
    if kk == base+nfpu+1,
        base = base+nfpu;
        count = count+1;
        nc = sprintf(['%0' ndigit 'd'],count);
        % load visual hull data
        save_name = [imgdir '/segment_part' nc '.mat'];
        if exist(save_name,'file')==2,
            load(save_name);
        else
            fprintf(1,'\nERROR: cannot find data ''%s''!\n',save_name);
            return;
        end;
    end;
    if active_images(kk),
        bk = kk-base;
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
            center = member_center(:,ii,kk);
            semiax = member_axes(:,ii,kk);
            vec = rodrigues(member_omega(:,ii,kk));
            ctt = center(:,ones(1,2));
            xyz = ctt+vec(:,1:2)*diag(semiax(1:2))*nlen_vec;
            for i=1:2,
                plot3([ctt(1,i);xyz(1,i)], [ctt(3,i);xyz(3,i)], -[ctt(2,i);xyz(2,i)],'color', palette(i,:), 'linewidth',2);
            end;
        end;
        title(['Segment result of frame: ', num2str(kk)]);
        axis equal vis3d off;
        view(az,el); hold off;
        pause(0.01);
        if check_direction,
            fprintf(1,'\nCheck axes direction of frame %d:\n',kk);
            flag = input('Some parts'' direciton need to be rectified? ([]=no, other=yes) ','s');
            if ~isempty(flag),
                flag = 1;
                while flag,
                    fprintf(1,'\nWhich parts need to rectify direciton?\nColor code for all members:\n');
                    for ii=1:nparts+1,
                        fprintf(1,['Part ' num2str(ii) ': ' color_cell{ii} ';\n']);
                    end;
                    idk = input(['Part number vector for rectification: (1~' num2str(nparts+1) ') ']);
                    if ~isempty(idk),
                        idk = round(idk);
                        flag = length(idk)>nparts+1 || any(idk<1) || any(idk>nparts+1);
                        % check indk if any color number is repeated
                        if ~flag,
                            for ii=1:nparts+1,
                                flag = sum(idk==ii)>1;
                                if flag, break; end;
                            end;
                        end;
                    end;
                    if flag,
                        fprintf(1,'\nUnexpected input! Please enter again!\n');
                    end;
                end;
                for ii= idk,
                    flag = 1;
                    while flag,
                        fprintf(1,'\nRectifying axes direciton of part %d (%s):\n',ii,color_cell{ii});
                        idk = input('Which axis is pointing inversely? (1=blue; 2=Green; 3=both;) ');
                        if ~isscalar(idk) || all(idk~=1:3),
                            fprintf(1,'Unexpected input! Please enter again!\n');
                            continue;
                        end;
                        flag = 0;
                    end;
                    temp = member_omega(:,ii,kk:end);
                    omega = reshape(temp,3,[]);
                    switch idk,
                        case 1,          % rodrigues(omega) times diag([-1,1,-1]) on the right hand side
                            temp(:) = trans_quat_axis(quatmul(trans_quat_axis(omega), [0; 1; 0; 0]));
                        case 2,          % rodrigues(omega) times diag([1,-1,-1]) on the right hand side
                            temp(:) = trans_quat_axis(quatmul(trans_quat_axis(omega), [1; 0; 0; 0]));
                        otherwise,      % rodrigues(omega) times diag([-1,-1,1]) on the right hand side
                            temp(:) = trans_quat_axis(quatmul(trans_quat_axis(omega), [0; 0; 1; 0]));
                    end;
                    member_omega(:,ii,kk:end) = temp;
                end;
            end;
        end;
    end;
end;

if check_direction,
    fprintf(1,'\nUpdating segment variables in ''segment_variables.mat''.\n');
    save_name = [imgdir '/segment_variables.mat'];
    string_save = ['save ' save_name ' n_frame ind_active active_images ndigit palette color_cell kparts nparts member_nax ' ...
        'member_center member_axes member_omega nfpu n_unit FramePS cyclePS FramePC Nstroke ind_reversal az el'];
    eval(string_save);
    fprintf(1,'done...\nNow you can refine the kinematics based on initial value and segmentation...\n');
end;

flag = input('Orientation of all members besides body look all right? ([]=no, other=yes) ','s');
if isempty(flag),
    fprintf(1,'\nRerun this script to correct directions of some members!\n');
    return;
end;
fprintf(1,'Now you can run script ''kinematics_extraction.m'' to extract kinematics ...\n');