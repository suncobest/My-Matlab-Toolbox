Instruction of the multiple camera 1D calibration toolbox

Requirements:
    MATLAB version: R2016b or later

Global variables definition:
    n_cam: number of cameras;
    n_ima: number of image frames;
    imsize: the 2*n_cam matrix showing the image size of all cameras;
    n_view = n_cam*n_ima;     % total number of images
    hand_list: the handedness of all cameras with respect to the reference coordinate system;
    hand_cc: the handedness of all cameras with respect to the main camera;
    np1D: number of points on the 1D calibration rod;
    npts = n_ima*np1D;        % total number of points
    rodlen: the 1D coordinates of np1D points on the rod;

    x_cell: the 2D pixel coordinates extracted from images stored in a vector cell, length of n_view,
            sequencing images after cameras. In each cell, np1D points are stored if the camera view have points,
            otherwise, the cell is left blank.
    Xrod: the 3D position of np1D points in all n_ima frames (ideal data);
    Xori: the origin of all rods in each frame (ideal data);
    Xdir: the direction of all rods in each frame (ideal data);
    Xp: the 3D position of np1D points in all n_ima frames (computed data);
    Xo: the origin of all rods in each frame (computed data);
    thph: the angular direction of all rods in each frame (theta and phi);

    active_imgviews: the n_cam*n_ima logical matrix showing which camera and image is active;
    active_images: the n_ima logical vector showing which frame have more than two cameras available images;
    ind_active: the index number of active images;

    A_cam: the n_cam*n_cam weighted adjacency matrix of all cameras; (smaller value means closer path)
    costs: the matrix of costs from each camera to other cameras;
    paths: the cell matrix to store the shortest path from each camera to other cameras;
    idm: the index number of the main camera;
    pathm: the cell vector which stores the shortest path from main camera to other cameras;

    est_fc_mat: switch to turn on/off the estimation of focal length (2*n_cam);
    center_optim_vec: switch to turn on/off the estimation of principal point (1*n_cam);
    est_dist_mat: switch to turn on/off the estimation of distortion coefficients (5*n_cam);
    est_alpha_vec: switch to turn on/off the estimation of pixel skew (1*n_cam);
    est_aspect_ratio_vec: switch to turn on/off the estimation of aspect ratio of focal length (1*n_cam);

    fc_mat: the 2*n_cam matrix to store the pixel focal length of all cameras;
    cc_mat: the 2*n_cam matrix to store the principle points of all cameras;
    alpha_vec: the 1*n_cam vector to store the skew coefficients of all cameras;
    kc_mat: the 5*n_cam matrix to store the distortion coefficients of all cameras;
    Omcw: the 3*n_cam matrix (rotation vector) to store the attitude of all cameras with respect to world reference frame;
    Tcw: the 3*n_cam matrix to store the position of all cameras with respect to world reference frame;
    Omcc: the 3*n_cam matrix (rotation vector) to store the attitude of all cameras with respect to the main camera;
    Tcc: the 3*n_cam matrix to store the position of all cameras with respect to the main camera;

    err_std0: the 2*1 reprojected error of all images before bundle adjustment (standard deviation);
    err_std: the 2*1 reprojected error of all images after bundle adjustment (standard deviation);
    err_cam: the 2*n_cam reprojected error of each camera after bundle adjustment (standard deviation);
    ex_max: the maximum reprojected error of all images after bundle adjustment;
    estdX: the reconstructed error of all points (standard deviation);
    errX_max: the maximum reconstructed error of points on the rod;

Main calibration subscripts:
    calib_1D_optim_func.m        % function version
    calib_1D_optim_multicam.m

Bundle adjustment function:
    optim_1D_multicam.m       % refine intrinsic and extrinsic parameters;
    optim_1D_extrinsic.m      % only refine extrinsic parameters;
    binocular_1D_optim.m      % pairwise intrinsic and extrinsic optimization;
    binocular_optim1D_extrinsic.m      % pairwise extrinsic optimization;

Input and output subroutine:
    simulate_1D_blocks.m      % generate 1D simulation data;
    write_simu_data.m         % convert MATLAB simulation data into "Goku" data structure;
    convert_goku_mat.m        % convert "Goku" system's output data format to my MATLAB data structure;
    plot_cameras_1D.m         % the drawing subscript showing all cameras and 3D points on rods;

Important tool function:
    comp_distortion.m         % compensate distortion with normalized 2D points; 
    apply_distortion.m        % apply distortion with normalized 2D points; 
    compose_motion2.m         % compose two Euclidian motion;
    compute_fundmatrix.m      % compute fundamental matrix using the 8-points algorithm;
    compute_rigid_refmotion.m       % compute the optimized rigid motion giving the feature points on a moving rigid (LM iteration);
    compute_rigid_refmotion2.m      % compute the optimized rigid motion giving the feature points on a moving rigid (R-T iteration);
    compute_Rt_pair.m         % compute the essential matrix using 8-points algorithm and decompose the essential matrix to R and T;
    project_points_mirror2.m  % projection of 3D world points to the 2D image plane;  
    compute_structure2.m      % reconstruct 3D points giving 2D correspondence and camera parameters;
    dijkstra.m                % the Dijkstra function to determine the shortest path giving the weighted adjacency matrix;
    gen_1D_points.m           % 3D points generator giving the rod origins and angular direction;
    distort_image.m           % add distortion to image giving the intrinsic parameters;
    rectify_image.m           % rectify distortion off an image giving the intrinsic parameters;
    normalise2dpts.m          % normalize 2D points with affine transformation (sqrt(2) of average side length);
    normalize2.m              % normalize pixel points to the unit focal length plane giving the intrinsic parameters;
    rigid_refmotion.m         % transform 3D coordinates under the Euclidian motion; 

Important math function:
    cartesian2spherical.m    % transform Cartesian coordinates into spherical coordinates;
    spherical2cartesian.m    % transform spherical coordinates into Cartesian coordinates;
    quatmean.m               % compute weighted quaternion means;
    rodrigues.m              % the Rodrigues function to convert rotation vector to orthogonal matrix and vice versa;
    dAB.m                    % compute the Jacobian matrix of product A*B with respect to A and B;
    trans_euler_mat.m        % transform Euler angle into rotation matrix;
    trans_quat_axis.m        % transform quaternion vectors into rotation vector;
    trans_quat_mat.m         % transform quaternion vectors into rotation matrix;

