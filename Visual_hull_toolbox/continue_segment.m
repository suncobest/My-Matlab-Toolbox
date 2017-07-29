plan = input('Which routine to follow? ([]=initialize and track, other=apply periodic constraint) ','s');
if isempty(plan),
    flag = 1;
    while flag,
        frame_list = input('Please input frame numbers for initialization and tracking: ');
        Nb = length(frame_list);
        if Nb==0,
            fprintf('Please enter again!\n')
            continue;
        elseif Nb>1,
            if ~isequal(diff(frame_list), ones(1,Nb-1)),
                fprintf('All frame numbers must be continuous! Please enter again!\n')
                continue;
            end;
        end;
        flag = 0;
    end;
    Not_ok = true(1,Nb);
    while any(Not_ok),
        % initiate and track the 1st motion cycle
        track_list = frame_list(Not_ok);
        track_members;
        check_segmentation;
        if any(Not_ok),
            disp('Some frames need to be rectified!');
            flag = input('Interpolate these frames or not? ([]=yes, other=no)');
            if isempty(flag),
                track_list = frame_list(Not_ok);
                ind = frame_list(~Not_ok);
                ind = ind(ind>track_list(1)-3 & ind<track_list(end)+3);
                Qt = NaN(4,nparts,length(track_list));
                for ii=1:nparts,
                    Qt(:,ii,:) = reshape(squad(ind, trans_quat_axis(reshape(member_omega(:,ii+1,ind),3,[])), track_list),4,1,[]);
                end;
                refine_body_wings;
                check_segmentation;
            end;
        end;
    end;
    
        Nstop = frame_list(end);
        for ii=1 : find(ind_maxspan<Nstop,1,'last'),
            jj = max(ind_maxspan(ii)-dfpc,1);
            ind = jj : min(ind_maxspan(ii)+dfpc, Nstop);
            vec = NaN(3,length(ind));
            for kk=ind,
                a = rodrigues(member_omega(:,2,kk));
                b = rodrigues(member_omega(:,3,kk));
                vec(:,kk-jj+1) = a(:,1)+b(:,1);
            end;
            [~,id] = min(sum(vec.^2,1));
            ind_maxspan(ii) = id+jj-1;
        end;
    
    for kk = frame_list,
        Q0 = trans_quat_axis(member_omega(:,:,kk));
        Q0(1:3,1) = -Q0(1:3,1);
        Qbw(:,:,kk) = quatmul(Q0(:,1), Q0(:,2:end));
    end;
    
else
    
    flag = 1;
    while flag,
        Nstop = input('Please enter the last frame number: ');
        if isempty(Nstop),
            fprintf('Please enter again!\n')
            continue;
        else
            Nb = find(ind_maxspan>Nstop,1);
            if Nb<3,
                fprintf('Unexpected input! Please enter again!\n')
                continue;
            end;
        end;
        flag = 0;
    end;
    fprintf(1,'\nThe function use orientation of wing wrt body in a cycle to estimate frames afterwards!\n');
    flag = input('Update reference period or not? ([]=no, other=yes) ','s');
    if ~isempty(flag),
        ind = round(linspace(ind_maxspan(Nb-2), ind_maxspan(Nb-1), max(8, round(FramePC/3))));
        t0 = (ind-ind(1))/(ind(end)-ind(1));
        Q0 = permute(Qbw(:,:,ind), [1,3,2]);
    end;
    
    period_sw = input('Apply periodic constraint to estimate wing''s: ([]=all direction, other=only incidence) ','s');
    period_sw = isempty(period_sw);
    flag = 1;
    while flag,
        func_flag = input('Which function to call: (1=ellipsoid_tracking, 2=ellipsoid_rotation, 3=ellipsoid_rotation2) ');
        if length(func_flag)~=1 || all(func_flag~=[1 2 3]),
            fprintf('Unexpected input! Please enter again!\n')
            continue;
        end;
        flag = 0;
    end;
    for nk = Nb:Nmspan,
        Nstart = Nstop+1;
        nstep = ind_maxspan(nk-1)-ind_maxspan(nk-2);
        idk = ind_maxspan(nk-1)+nstep;
        if nk==Nmspan,
            Nstop = n_frame;
        else
            Nstop = idk+dfpc;
        end;
        track_list = Nstart:Nstop;
        tn = mod(track_list-ind_maxspan(nk-1),nstep)/nstep;
        for ii=1:nparts,
            Qbw(:,ii,track_list) = reshape(squad(t0, Q0(:,:,ii), tn),4,1,[]);
        end;
        % estimate direction by the last period
        track_body_wings;
        % update image number of maximum span
        jj = idk-dfpc;
        ind = jj : min(idk+dfpc, Nstop);
        vec = NaN(3,length(ind));
        for ii=ind,
            a = rodrigues(member_omega(:,2,ii));
            b = rodrigues(member_omega(:,3,ii));
            vec(:,ii-jj+1) = a(:,1)+b(:,1);
        end;
        [~,id] = min(sum(vec.^2,1));
        ind_maxspan(nk) = id+jj-1;
    end;
end;

fprintf(1,'\nSaving segment variables in ''segment_variables.mat''.\n');
save_name = [imgdir '/segment_variables.mat'];
string_save = ['save ' save_name ' n_frame ind_active active_images ndigit palette color_cell kparts nparts '...
    'member_center member_axes member_omega nfpu n_unit FramePS cyclePS FramePC Nstroke nfpc ind_down '...
    'ind_up ind_maxdown ind_maxup ind_maxspan wroot Qbw dfpc az el'];
eval(string_save);
fprintf(1,'done...\n');
