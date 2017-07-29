I=uint8(zeros(300,400*3,3));
Backpath = 'H:\昆虫飞行\三机标定拍摄\14.01.15\selected\pronation and supination\front\processed\';
Rightpath = 'H:\昆虫飞行\三机标定拍摄\14.01.15\selected\pronation and supination\right\processed\';
Toppath = 'H:\昆虫飞行\三机标定拍摄\14.01.15\selected\pronation and supination\top\processed\';
fmt = 'png';
savepath = 'H:\昆虫飞行\三机标定拍摄\14.01.15\selected\pronation and supination\3view movie\';


Backdir = dir([Backpath '*.' fmt]);
Rightdir = dir([Rightpath '*.' fmt]);
Topdir = dir([Toppath '*.' fmt]);

nframe = length(Backdir);
for i = 1:nframe,
    Back=imread([Backpath Backdir(i).name]);
    Right=imread([Rightpath Rightdir(i).name]);
    Top=imread([Toppath Topdir(i).name]);
    I(:,1:400,:)=Back;
    I(:,401:800,:)=Right;
    I(:,801:end,:)=Top;
    % image(I)
    % axis image
    imname = sprintf('3view_%03d',i);
    imwrite(I,[savepath imname '.png'],'png');
end;
    

   