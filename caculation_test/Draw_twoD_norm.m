function Draw_twoD_norm
I = zeros(100);
mu = [20,30];
P = [100,0;0,100];
for i = 1:100
    for j=1:100
        I(i,j)=1/sqrt((2*pi)^2*det(P))*exp(-[i-mu(1),j-mu(2)]/P*[i-mu(1),j-mu(2)]'/2);
    end
end
I = I/I(mu(1),mu(2));

set(0,'DefaultFigureWindowStyle','docked') 
imshow(I, 'InitialMagnification',500, 'Border','tight')





% % new figure
% hFig = figure;
% 
% % try show image at full size (suppress possible warning)
% s = warning('off', 'Images:initSize:adjustingMag');
% imshow(I, 'InitialMagnification',500, 'Border','tight')
% warning(s);
% 
% % handle figure resize events
% hAx = gca;
% set(hFig, 'ResizeFcn',{@onResize,hAx})
% 
% % call it at least once
% feval(@onResize,hFig,[],hAx);
% 
% % enable panning tool
% pan on
% end
% 
% function onResize(o,e,hAx)
% % get axes limits in pixels
% oldUnits = get(hAx, 'Units');    % backup normalized units
% set(hAx, 'Units','pixels')
% pos = get(hAx, 'Position');
% set(hAx, 'Units',oldUnits)       % restore units (so it auto-resize)
% 
% % display the top left part of the image at magnification 100%
% xlim(hAx, [0 pos(3)]+0.5)
% ylim(hAx, [0 pos(4)]+0.5)
% end