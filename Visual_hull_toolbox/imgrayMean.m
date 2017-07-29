function im = imgrayMean(dirname,basename,fmt)

% dirname is the directory of images to be processed
% basename is the base of image name without number and format
% the class of IM is double.

imfile = dir([dirname '/' basename '*.' fmt]);  % list of image files
nima = length(imfile);  % number of basename*.fmt
if nima == 0,
    im = [];
    fprintf(1,'No specified image found in the directory!\n');
    return;
end;

im = double(imread([dirname '/' imfile(1).name])); % initialization
np = size(im,3);
if np==3,
    im = 0.299*im(:,:,1) + 0.587*im(:,:,2) + 0.114*im(:,:,3);
end;
if nima>1,
    for kk = 2:nima,
        I = double(imread([dirname '/' imfile(kk).name]));
        if np==3,
            I = 0.299* I(:,:,1) + 0.587*I(:,:,2) + 0.114*I(:,:,3);
        end;
        im = im + I;
    end;
    im = im/nima;
end;

% figure;
% image(uint8(im));
% colormap(gray(256));
% axis image;
return;
