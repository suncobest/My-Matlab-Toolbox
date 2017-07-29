I=imread('D:\My Documents\MATLAB\visionbook\01Intro\images\Ondra_sampling.jpg');

map = gray(256);

minI = min(I(:));
maxI = max(I(:));

Id = 255*(I - minI)/(maxI - minI);

figure(2);
image(Id);
colormap(map);


for q=0:20:100
    filename=sprintf('quality_%03d.jpg',q);
    imwrite(I, filename,'quality',q);
end

%%
x = 0:0.01:1;
figure
filename = 'testAnimated.gif';

for n = 1:0.5:5
y = x.^n;
plot(x,y)
drawnow
frame = getframe(1);
im = frame2im(frame);
[A,map] = rgb2ind(im,256); 
	if n == 1;
		imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
	else
		imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
	end
end


% eval(['y = y_' num2str(kk) ';']);
% xc = cornerfinder(y+1,I,winty,wintx); % the four corners
% eval(['wintx_' num2str(kk) ' = wintx;']);
% eval(['winty_' num2str(kk) ' = winty;']);
% eval(['x_' num2str(kk) '= xc - 1;']);

% title(['Image ' num2str(kk) ' - Image points (+) and reprojected grid points (o)']);
% eval(['plot(x_' num2str(kk) '(1,:)+1,x_' num2str(kk) '(2,:)+1,''r+'');']);
% eval(['plot(y_' num2str(kk) '(1,:)+1,y_' num2str(kk) '(2,:)+1,''' colors(rem(kk-1,6)+1) 'o'');']);
% eval(['quiver(y_' num2str(kk) '(1,:)+1,y_' num2str(kk) '(2,:)+1,ex_' num2str(kk) '(1,:),ex_' num2str(kk) '(2,:),1,''' colors(rem(kk-1,6)+1) ''');']);
