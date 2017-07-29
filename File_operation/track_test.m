% Locate the image
[imname,pathname]= uigetfile({'*.jpg;*.JPG;*.bmp;*.BMP;*.png;*.PNG;*.gif;*.GIF;*.tif','Image Files (*.jpg,*.bmp,*.png,*.gif,*.tif)';...
    '*.jpg','JPEG Files(*.jpg)'; '*.bmp','BMP Files(*.bmp)'; '*.png','PNG Files(*.png)'; '*.gif','GIF Files(*.gif)';...
    '*.tif','TIFF Files(*.tif)'; '*.*','All Files(*.*)'},'Select the first frame image!');    % select path for image
if imname == 0
    fprintf(1,'No image selected!\n\n');
    return;
end;

I = imread(fullfile(pathname,imname));
if ndims(I)==3,
    I =  0.299 * I(:,:,1) + 0.587 * I(:,:,2) + 0.114 * I(:,:,3);
end;
I0 = I;

imdir = dir(fullfile(pathname,['*.' imname(end-2:end)]));
nframe = length(imdir);
ndigit = floor(log10(nframe));
start_numstr = imname(end-4-ndigit:end-4); 


imname = imdir(1).name;
start = str2double(start_numstr) - str2double(imname(end-4-ndigit:end-4))+1;
if exist('start_num','var'),             % in case that the images corners have already been extracted with different first frame
    start_idx = start - start_num + 1;
end;
start_num = start;

figure(1);
image(I);
colormap(gray(256));
set(1,'color',[1 1 1]);
axis image;

% extract corner points
if ~exist('npoints','var'),
    npoints = input('How many corner points do you want to extract? ( [] = 1 )');
    if isempty(npoints),
        npoints = 1;
    end;
end;

% set the size of window for finding corners
wint_default = 7;
disp('Window size for corner finder (wintx and winty):');
if ~exist('wintx','var') || ~exist('winty','var'),
    wintx = input(['wintx ([] = ' num2str(wint_default) ') = ']);
    if isempty(wintx), wintx = wint_default; end;
    wintx = round(wintx);
    winty = input(['winty ([] = ' num2str(wint_default) ') = ']);
    if isempty(winty), winty = wint_default; end;
    winty = round(winty);
end;
fprintf(1,'Window size = %d x %d\n',2*wintx+1,2*winty+1);

figure(1);
image(I);
colormap(gray(256));
set(1,'color',[1 1 1]);
axis image;
hold on;

if ~exist('start_idx','var') || start_idx < 1,
    title(['Click ' num2str(npoints) ' points on the image: ' start_numstr]);    
    x= [];y = [];
    for count = 1:npoints,
        [xi,yi] = ginput(1);
        xxi = cornerfinder([xi;yi],I,winty,wintx);
        xi = xxi(1);
        yi = xxi(2);
        x = [x;xi];
        y = [y;yi];
        figure(1);
        plot(xi,yi,'g+','linewidth',2,'Markersize',10);
        plot(xi + [wintx+.5, -(wintx+.5), -(wintx+.5), wintx+.5, wintx+.5],yi + [winty+.5, winty+.5, -(winty+.5), -(winty+.5), winty+.5],'-','color',[ 1.000 0.314 0.510 ],'linewidth',2);
        drawnow;
    end;  
    xt = [x,y]';
else
    xt = xtrack(start_idx*2-1:start_idx*2,:);
    figure(1);
    title(['Corner points on the image: ' start_numstr]);  
    plot(xt(1,:),xt(2,:),'r+','linewidth',2,'Markersize',10);
end;
hold off;


% tracking process
nn = 3;
wxtrack = wintx*nn;
wytrack = winty*nn;

inspeed = zeros(size(xt)); 
xtrack = xt;

for i = start_num+1:nframe,
    imname = imdir(i).name;
    I = imread(fullfile(pathname,imname));
    if ndims(I)==3,
        I =  0.299 * I(:,:,1) + 0.587 * I(:,:,2) + 0.114 * I(:,:,3);
    end;
    
    [xc,good,outspeed] = cornertracker(xt,I0,I,wintx,winty,wxtrack,wytrack,inspeed);
 
    % show image in the first window
    figure(1);
    image(I);
    colormap(gray(256));
    set(1,'color',[1 1 1]);
    axis image;
    hold on;
    title(['Corner points on the image: ' imname(end-4-ndigit:end-4)]);  
    plot(xc(1,:),xc(2,:),'g+','linewidth',2,'Markersize',10);
    hold off;
    
    I0 = I;
    xt = xc;
    xtrack = [xtrack; xt];
    inspeed = outspeed;
    
%     pause(0.02);
end;

trackdata = fullfile(pathname,'xtrack.txt');
if exist(trackdata,'file') ~= 2,
    fid=fopen(trackdata,'w+');  % 'r+','w+','a+'
    [xm,xn]=size(xtrack);
    for i = 1:xm,
        for j = 1:xn,
            fprintf(fid,'%10.4f\t',xtrack(i,j));
        end;
        fprintf(fid,'\r\n');
    end
    fclose(fid);
end;

% xtrack = load(trackdata);  % read tracking data