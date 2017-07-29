% Locate the image
[imname,pathname]= uigetfile({'*.jpg;*.JPG;*.bmp;*.BMP;*.png;*.PNG;*.gif;*.GIF;*.tif','Image Files (*.jpg,*.bmp,*.png,*.gif,*.tif)';...
    '*.jpg','JPEG Files(*.jpg)'; '*.bmp','BMP Files(*.bmp)'; '*.png','PNG Files(*.png)'; '*.gif','GIF Files(*.gif)';...
    '*.tif','TIFF Files(*.tif)'; '*.*','All Files(*.*)'},'Select the first image for tracking!');    % select path for image
if imname == 0
    fprintf(1,'No image selected!\n\n');
    return;
end;

trackdata = fullfile(pathname,'xtrack.txt');
if exist(trackdata,'file')~=2,
    fprintf(1,'No tracking data (xtrack.txt) was founnd in the directory!\n\n');
    return;
end;
xtrack = load(trackdata);  % read tracking data

imdir = dir(fullfile(pathname,['*.' imname(end-2:end)]));
nframe = length(imdir);
ndigit = floor(log10(nframe));
start_numstr = imname(end-4-ndigit:end-4); 

imname = imdir(1).name;
start_num = str2double(start_numstr) - str2double(imname(end-4-ndigit:end-4))+1;

%%

for i = start_num : nframe,
    imname = imdir(i).name;
    xc = xtrack(i*2-1:i*2,:);          % the corner points for the i-th frame
    
    I = imread(fullfile(pathname,imname));
    if ndims(I)==3,
        I =  0.299 * I(:,:,1) + 0.587 * I(:,:,2) + 0.114 * I(:,:,3);
    end;
    
    % show image in the second window
    figure(1);
    image(I);
    colormap(gray(256));
    set(1,'color',[1 1 1]);
    axis image;
    hold on;
    title(['Corner points on the image: ' imname(end-4-ndigit:end-4)]);
    plot(xc(1,:),xc(2,:),'g+','linewidth',2,'Markersize',10);
    hold off;
    
    pause(0.1);
end;