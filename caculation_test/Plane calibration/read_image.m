% Locate the image you want to calibrate
if exist('imname','var') && exist('I','var') && ismatrix(I)
    fprintf(1,'Do you want to reload image? \n')
    reply = input('yes or no? ([] = yes) >>','s');
    if isempty(reply)
        reply = 'yes'; 
    end
    if strcmp(reply, 'no')
        fprintf(1,'Program load the image that already exist!\n')
        figure(2);
        image(I);
        colormap(gray(256));
        set(2,'color',[1 1 1]);
        axis image;
        return
    elseif ~strcmp(reply,'yes')
        error('Unexpected input!')
    end
end

[imname,pathname]= uigetfile({'*.jpg;*.bmp;*.png;*.gif;*.tif','Image Files (*.jpg,*.bmp,*.png,*.gif,*.tif)';...
    '*.jpg','JPEG Files(*.jpg)'; '*.bmp','BMP Files(*.bmp)'; '*.png','PNG Files(*.png)'; '*.gif','GIF Files(*.gif)';...
    '*.tif','TIFF Files(*.tif)'; '*.*','All Files(*.*)'},'Select the first frame image!');    % select path for image
if imname == 0
    fprintf(1,'No image selected!\n\n')
    return
end

% Read image from pathname to the matrix I
I = imread(fullfile(pathname,imname));
fprintf(1,'Image loaded!\n\n');


% Convert RGB to Gray
if ndims(I)==3
    I =  0.299 * I(:,:,1) + 0.587 * I(:,:,2) + 0.114 * I(:,:,3);
end
% I_size = size(I);

[ny,nx] = size(I);
% show image in the second window
figure(2);
image(I);
colormap(gray(256));
set(2,'color',[1 1 1]);
axis image;


