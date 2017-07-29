function edge = invprj_vector(pedge,dir_plane)
% pedge is the vector in pixel
% dir_plane is the direction of the calibrated plane
% edge is the physical vector corresponding to pedge in the scene (in pixel)
dir_prj = dir_plane(1:2);
dir_prj = [dir_prj/norm(dir_prj),0];
pedge_dir = sum(pedge.*dir_prj);
edge = pedge_dir;
pedge_dir = pedge_dir * dir_prj;
pedge_ndir = pedge - pedge_dir;
ndir_prj = pedge_ndir/norm(pedge_ndir);   % base vectors (dir_prj, ndir_prj, [0 0 1]) in the picture 
hand = cross(dir_prj, ndir_prj);          % with the counterpart base vectors (hand*cross(ndir_prj,dir_plane), ndir_prj, dir_plane) in the scene
hand = hand(3);
edge = (edge/dir_plane(3)* hand * cross(ndir_prj,dir_plane) + pedge_ndir);

return