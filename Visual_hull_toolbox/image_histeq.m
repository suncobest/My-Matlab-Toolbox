function [out,list] = image_histeq(img)
% IMAGE_HISTEQ Enhance contrast using histogram equalization.
%   IMAGE_HISTEQ enhances the contrast of images by transforming the values in an
%   intensity image. out is double, need to transform to uint8 before saving. See also histeq.

if islogical(img),
    out = img;
    fprintf(1,'Logical picture detected! No need to run this funciton!\n');
    return;
end;

[ny,nx,nc] = size(img);
img = double(img);
if nc==3,
    img = 0.299 *img(:, :,1)+ 0.587*img(:,:,2) + 0.114*img(:,:,3);
end;

nbins = 0:255;
[nn, ind] = histc(img(:),nbins);
list = cumsum(nn);
list = round(list*255/(nx*ny));             % new gray scale map
out = reshape(list(ind),ny,nx);

return;



%% Test
I = imread('tire.tif');
J = uint8(image_histeq(I));
figure, imshowpair(I,J,'montage')