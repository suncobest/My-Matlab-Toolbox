% 3D data are transformed into specific system: X1 = R1*H1*Xw+T1;
om1 = Omcw(:,1);
R1t = rodrigues(-om1);
T1 = Tcw(:,1);
hand1 = hand_list(:,1);
handcc = hand1*hand_list;
XX = rigid_refmotion(Xrod,om1,T1,hand1);

save_name = 'X3d';
fid = fopen([save_name '.txt'],'wt');
fprintf(1,['\nSaving 3D data in ''' save_name '.txt''\n']);
n0 = zeros(1,np1D*3);
active_view = any(active_imgviews,1);
for kk=1:n_ima,
    ii = (kk-1)*np1D;
    if active_view(kk),
        x_kk = XX(:,ii+1:ii+np1D);
        fprintf(fid,'%d',kk-1);
        fprintf(fid,' %.6f',x_kk(:));
        fprintf(fid,'\n');
    else
        fprintf(fid,'%d',kk-1);
        fprintf(fid,' %d',n0);
        fprintf(fid,'\n');
   end;
end;
fclose(fid);

n0 = zeros(1,np1D*2);
n1 = -ones(1,np1D*2);
save_name = 'camcfg';
fid = fopen([save_name '.ini'],'wt');
fprintf(1,['\nSaving camera parameters in ''' save_name '.ini''\n']);
active_view = sum(active_imgviews,2)>10;
for pp=1:n_cam,
    if active_view(pp);
        fc = fc_mat(:,pp);
        cc = cc_mat(:,pp);
        alpha_c = alpha_vec(pp);
        kc = kc_mat(:,pp);
        handkk = handcc(pp);
        % Calibration matrix:
        KK = [fc(1) fc(1)*alpha_c cc(1);0 fc(2) cc(2); 0 0 1];
        omwkk = Omcw(:, pp);
        Twkk = Tcw(:, pp);
        [omck, Tck] = compose_motion2(-om1,-R1t*T1,omwkk,Twkk,handkk);  %  inverse composition
        % save cameras
        fprintf(fid,'[camera_%03d]\nKK =',pp-1);
        fprintf(fid,' %.8f',KK(:));
        fprintf(fid,'\nkc =');
        fprintf(fid,' %.8f',kc(:));
        fprintf(fid,'\nimsize =');
        fprintf(fid,' %d',imsize(:,pp));
        fprintf(fid,'\nOmc =');
        fprintf(fid,' %.8f',omck);
        fprintf(fid,'\nTc =');
        fprintf(fid,' %.8f',Tck);
        fprintf(fid,'\nhandedness = %d\n\n',handkk);

        % save 2D data
        save_name = sprintf('x2d_cam%03d',pp-1);
        ind = fopen([save_name '.txt'],'wt');
        fprintf(1,['\nSaving 2d pixel data of camera %d in ''' save_name '.txt''\n'],pp);
        save_name = sprintf('x2dn_cam%03d',pp-1);
        id = fopen([save_name '.txt'],'wt');
        fprintf(1,['\nSaving 2d normalized data of camera %d in ''' save_name '.txt''\n'],pp);
        for kk=1:n_ima,
            if active_imgviews(pp,kk),
                kth = (kk-1)*n_cam+pp;
                x_kk = x_cell{kth};
                xn = normalize_pixel(x_kk,fc,cc,kc,alpha_c);
                fprintf(ind,'%d',kk-1);
                fprintf(ind,' %.6f',x_kk(:));
                fprintf(ind,'\n');
                fprintf(id,'%d',kk-1);
                fprintf(id,' %.6f',xn(:));
                fprintf(id,'\n');
            else
                fprintf(ind,'%d',kk-1);
                fprintf(ind,' %d',n1);
                fprintf(ind,'\n');
                fprintf(id,'%d',kk-1);
                fprintf(id,' %d',n0);
                fprintf(id,'\n');
            end;
        end;
        fclose(ind);
        fclose(id);
    end;
end;

fclose(fid);
