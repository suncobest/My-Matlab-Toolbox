% swap x and y coordinates of corner points
% This script apply to transparent checkerboard
% n_sq_x与n_sq_y表示四个顶点内x和y方向上的矩形数量

if ~exist('x_cell','var'),
    fprintf(1,['\nThere is no corner points data to process!\n'...
        'Please load ''multimirror_calib_data.mat'' or ''multicam_calib_data.mat"...\n']);
    return;
end;
pp = input(['Which view need to swap x-y corrdinates? ([' num2str(1:n_cam) '])']);
assert(length(pp)==1 && any(pp==(1:n_cam)),'Unexpected input! Please choose one view at a time!');
imas = input(['Range of image number to process: (1~' num2str(n_ima) ')']);
assert(isvector(imas) && min(imas)>=1 && max(imas)<=n_ima,'Unexpected input for image numbers!');

for kk=imas,
    kth = (kk-1)*n_cam+pp;
    x = x_cell{kth};
    if ~isempty(x),
        n_sq_x = n_sq_mat(1,kth);
        n_sq_y = n_sq_mat(2,kth);
        Np = (n_sq_x+1)*(n_sq_y+1);
        % generate original array index, four corner (y,xy,x,o)
        % 排成点阵,行向为y方向，列向为x方向，第一个点在右上角，沿x列向排列，原点在左上角
        idx=fliplr(reshape(1:Np, n_sq_x+1, n_sq_y+1));

        % transform corner array index
        idx=rot90(idx,-1);   % 将点阵顺时针旋转90度，将最底下一行变成第一列
        idx=idx(:)';  % 得到变换后的点阵索引

         % final indexing
        x_cell{kth}=x(:,idx);

        % check reproject_corners_mirror or reproject_corners_multicam for
        % definition of axes X and Y.
        dX= dXY_mat(1, kth) ;
        dY = dXY_mat(2, kth);
        n_sq_mat(1,kth)= n_sq_y;
        n_sq_mat(2,kth)= n_sq_x;
        Xi = reshape(((0:n_sq_y)*dX)'*ones(1,n_sq_x+1),1,Np);
        Yi = reshape(ones(n_sq_y+1,1)*(n_sq_x:-1:0)*dY,1,Np);
        X_cell{kth} = [Xi;Yi];
    end;
end;

fprintf(1,['\nPlease save rectified data to ''multimirror_calib_data.mat'' or ''multicam_calib_data.mat" if all error is done!\n']);
