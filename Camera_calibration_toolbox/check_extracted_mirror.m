for kk =  1:n_ima,
    for pp = 1:n_cam,
        kth = (kk-1)*n_cam+pp;
        x = x_cell{kth};
        if isempty(x) || isnan(x(1)),   % 若x为[]或NaN，则此视角无效；
            active_imgviews(pp,kk) = 0;
            x_cell{kth} = [];
            X_cell{kth} = [];
            dXY_mat(:, kth) = NaN(2,1);
            n_sq_mat(:, kth) = NaN(2,1);
        end;
    end;
    if all(active_imgviews(:,kk)==0),
        fprintf(1,['WARNING: Not a single view have grid corners on image ' ...
                   num2str(kk) ' - This image is now set inactive!\n']);
    end;
end;
active_images = any(active_imgviews,1);
ind_active = find(active_images);
