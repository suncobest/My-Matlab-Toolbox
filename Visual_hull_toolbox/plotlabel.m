function plotlabel(frame,n)
% Using the function BWLABEL to label different objects: [ob_label,n_ob] = bwlabel(frame,n)
% frame is logical image; n is connectivity; ob_label is the label matrix
% the same size of frame; 0 means background; n_ob is number of objects identified. 
% plot BW images with different colors or numbers according to objects' label.

assert(islogical(frame) && ismatrix(frame),'The first argument must be logical matrix!');
if nargin<2,
    n=8;
else 
    assert(~isempty(n) && (n==4 || n==8),'The2nd argument (connectivity) must be 4 or 8!');
end;

if issparse(frame),
    frame = full(frame);
end;

[ob_label,n_ob] = bwlabel(frame,n);  % 4 or 8-connected (n=4 ro 8)
colors = lines(n_ob);   % lines is a color space function

figure(1);
imR = ones(size(frame));
imG = imR;
imB = imR;
for label = 1:n_ob,
    logicId = (ob_label==label);  % logical indexing
    imR(logicId) = colors(label,1);  % fill object with colors
    imG(logicId) = colors(label,2);
    imB(logicId) = colors(label,3);
end;

im = cat(3,imR,imG,imB);  % Concatenate arrays
image(im2uint8(im));   % image函数的参数类型最好为uint8,16
set(1,'color',[1 1 1]);
hold on;
for label = 1:n_ob,
    [yi,xi] = find(ob_label==label);
    text(mean(xi),mean(yi),num2str(label),'color','k');
end
hold off;

return;




colors = 'bgrcmy';

if  nxy(2) > 100 || nxy(1) > 80
    imR = ones(nxy);
    imG = imR;
    imB = imR;
    for label = 1:n_ob
        logicId = (ob_label==label);  % logical indexing
        imR(logicId) = colors(label,1);  % fill object with colors
        imG(logicId) = colors(label,2);
        imB(logicId) = colors(label,3);    
    end
    im = cat(3,imR,imG,imB);  % Concatenate arrays
    image(im2uint8(im));
    set(1,'color',[1 1 1]);   
else
    image(im2uint8(~frame));  % image函数的参数类型最好为uint8,16
    colormap(gray(256));
    set(1,'color',[1 1 1]);
    hold on;
    for label = 1:n_ob
        [yi,xi] = find(ob_label==label);
        text(xi,yi,num2str(label),'color',colors(label,:));
        %         plot(xi,yi,'.','color',colors(rem(label-1,6)+1));
    end
end

hold off;