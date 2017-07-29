% This script make gif of image sequence
% initialization of number of image frames
n_frame = 0;
if exist('imgdir','var') && ~isempty(imgdir),
  fprintf(1,'Run the script ''imgdir = [];'' to redirect image sequence path!\n');
  if ~exist('imgbase','var') || isempty(imgbase),
      imgbase = input('Basename of image frames (without number nor suffix): ','s');
  else
      fprintf(1,'Run the script ''imgbase = [];'' to change image basename!\n');
  end;
  if ~exist('imgfmt','var') || isempty(imgfmt),
      imgfmt = input('Format of image frames (only suffix): ','s');
  else
      fprintf(1,'Run the script ''imgfmt = [];'' to change image format!\n');
  end;
  %file name prefix number of original image
  imPreNum = [imgdir '/' imgbase];
  [n_frame, strnum_frame] = check_image_sequence(imPreNum, imgfmt);
end;

if n_frame ==0,
  fprintf(1,'No images found! Please relocate!\n');
end;

while n_frame == 0,
[imgbase,imgdir]= uigetfile({'*.jpg;*.bmp;*.png;*.gif;*.tif','Image Files (*.jpg,*.bmp,*.png,*.gif,*.tif)';
                              '*.jpg','JPEG Files(*.jpg)';'*.bmp','BMP Files(*.bmp)'; '*.png','PNG Files(*.png)';
                              '*.gif','GIF Files(*.gif)'; '*.tif','TIFF Files(*.tif)';
                              '*.*','All Files(*.*)'},'Select the first frame image!');    % 选择图片路径
if imgbase == 0,
  n_frame = 0;
  fprintf(1,'Please relocate images!\n');
else
  loc_ext = strfind(imgbase,'.');
  imgfmt = imgbase(loc_ext(end)+1:end);   % format of images
  imgbase = input('Basename of image frames (without number nor suffix): ','s');
  if imgdir(end)=='\' || imgdir(end)=='/',
    imgdir(end) = [];
  end;
  imPreNum = [imgdir '/' imgbase];
  [n_frame, strnum_frame] = check_image_sequence(imPreNum, imgfmt);
end;
end;

save_name = [imgdir '/play_' imgbase '.gif'];
if exist(save_name,'file')==2,
  fprintf(1,'\nGIF file ''play_%s.gif'' detected in the directory ...\n',imgbase);
  flag = input('Do you want to overwrite it or not? ([]=no, other=yes) ','s');
  if isempty(flag),
      fprintf(1,'\nGIF saving process have been terminated ...\n');
      return;
  else
      delete(save_name);
  end;
end;

if n_frame>80, 
  step = ceil(n_frame/50);
  fprintf(1,'There are too many images to make a gif!\nSetting skip ratio to %d ...\n',step);
else
  step = 1;
end;

wt = 0.1;
height = 3;   % unit: inch
% string form of image number digit
ndigit = num2str(floor(log10(n_frame))+1);
frame_kk = imread([imPreNum strnum_frame{1} '.' imgfmt]);
[ny,nx,nc] = size(frame_kk);
tpy = ny/20;
tpx = nx-tpy;  % nx/2;
resolution = round(ny/height);   % dpi

% [x,y,w,h], i.e. left top corner and width and height
mask = [nx*2/3, 0, nx/3, ny/8];
improc = 24:61; %1:step:n_frame,
for kk = improc,
  frame_kk = imread([imPreNum strnum_frame{kk} '.' imgfmt]);
  figure(2);
  image(frame_kk);
  if nc==1,
      if islogical(frame_kk),
          colormap(gray(2));
      else
          colormap(gray(256));
      end;
  end;
  hold on;

  % mask to cover old text
  % rectangle('Position',mask,'facecolor','k','edgecolor','none','linewidth',1);
  % framenb = sprintf(['%0' ndigit 'd'],kk-improc(1)+1);
  % text(tpx,tpy,['\it\fontname{Arial}\fontsize{16}\color{white}Frame = ' ...
  %     '\color{yellow}',framenb],'HorizontalAlignment','right');       %  'center');

%         text(tpx,tpy,['\it\fontname{Arial}\fontsize{16}\color{black}Frame = ' ...
%         '\color{red}',framenb],'HorizontalAlignment','right');
  axis off equal tight;
  % 设置绘图区域（坐标系）在窗口中的位置, position = [left, bottom, width, height]
  set(gca,'position',[0 0 1 1]);
%     set(gcf, 'unit', 'normalized', 'position', [0.1,0.1,0.8,0.8]);  % 设置窗口全屏
    set(gcf,'PaperPositionMode','Auto','PaperUnits','inches','PaperPosition',[0 0 nx ny]/resolution);
    hold off;
    drawnow;
    frame = getframe(gcf);
    im = frame2im(frame);
    [I,map] = rgb2ind(im,256);
    
    if  kk == improc(1),
        imwrite(I,map,save_name,'gif','LoopCount',Inf,'DelayTime',wt);
    else
        imwrite(I,map,save_name,'gif','WriteMode','append','DelayTime',wt);
    end;
end;

fprintf(1,'\nGIF file has been saved as ''play_%s.gif'' in the image directory!\n',imgbase);
