%% add time to pics
% 1 设置paperposition为manual
% set(gcf,'PaperPositionMode', 'manual')
% [ auto | {manual} ]
% 
% 2 设置paperunit   
% set(gcf,'PaperUnits','inches')
% [ {inches} | centimeters | normalized | points ]
% 
% 3 设置paperposition
% set(gcf,'PaperPosition',[left,bottom,width,height])

% % default fontsize 12


width = 800;      % unit: points
height = 700;
R = height/2;       % resolution: dpi
tpx = width/2;    % unit: points
tpy = 50;

I=im2uint8(zeros(height,width,3));
fpath = 'L:\InsectFlying\Mayun\MM\3';
fmt = 'png';
savepath = 'L:\InsectFlying\Mayun\movie3';

picdir = dir(fullfile(fpath,['*.' fmt]));
I0 = imread(fullfile(fpath,picdir(1).name));
[ny,nx,~] = size(I0);

nframe = length(picdir);     % 1;  
for i = 1:nframe,
    figure(1);
    image(I);
    axis off image;
    hold on;
    strnb = sprintf('%02d',i);
    save_name = fullfile(savepath,[strnb '.' fmt]);
    text(tpx,tpy,['\it\fontname{Arial}\fontsize{10.5}\color{white}Frame = \color{blue}', strnb],'HorizontalAlignment','center');
    set(gca,'position',[0 0 1 1]); % 设置绘图区域（坐标系）在窗口中的位置, position = [left, bottom, width, height]
    set(gcf,'PaperPositionMode','manual','PaperUnits','inches','PaperPosition',[0 0 width height]/R);
    drawnow;
    print(['-d' fmt],['-r' num2str(R)],save_name);          % print('-dpng','-r200',fnames(i))
    hold off;
    
%     keyboard;           % switch to turn on/off image adding
    I0= imread(fullfile(fpath,picdir(i).name));
    It = imread(save_name);
    It(end-ny+1:end,:,:)=I0;
    imwrite(It,save_name,fmt);
end;



%% add text and time to pics
% see connect_views.m

I=im2uint8(ones(400,400*3,3));
fpath = 'G:\Test Data\14.01.15\selected\pronation and supination\3view movie\images\original';
fmt = 'png';
savepath = 'G:\Test Data\14.01.15\selected\pronation and supination\3view movie\images\marked';

picdir = dir(fullfile(fpath,['*.' fmt]));

nframe = length(picdir);
for i = 1:nframe,
    I0=imread(fullfile(fpath,picdir(i).name));
    I(101:end,:,:)=I0;
    
    figure(1);
    image(I);
    axis off image;
    hold on;
    strnb = sprintf('%02d',i);
    text(600,25,['\it\fontname{Arial}\fontsize{24}\color{black}frame = \color{red}', strnb],'HorizontalAlignment','center');
    text(200,75,'\bf\fontname{Arial}\fontsize{24}\color{black}Back','HorizontalAlignment','center');
    text(600,75,'\bf\fontname{Arial}\fontsize{24}\color{black}Right','HorizontalAlignment','center');
    text(1000,75,'\bf\fontname{Arial}\fontsize{24}\color{black}Top','HorizontalAlignment','center');
    set(gcf, 'unit', 'normalized', 'position', [0,0.1,1,0.535]);  % 设置窗口大小，宽全屏，高为0.435倍的屏高
    set(gca,'position',[0 0 1 1]); % 设置绘图区域（坐标系）在窗口中的位置, position = [left, bottom, width, height]
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    
    imwrite(A,map,fullfile(savepath,[strnb '.png']),'png');
    hold off;
   
end;
    