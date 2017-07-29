function [intensityC,label]=kmeans_segment(ima,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   kmeans image segmentation
%
%   Input:
%          ima: grey color image
%          k: Number of classes
%   Output:
%          intensityC: vector of intensity center 
%          label: clasification image label
%
%   Author: Jose Vicente Manjon Herrera
%    Email: jmanjon@fis.upv.es
%     Date: 27-08-2005
%  Modified by Pengfei Zhang
%  Date: 27-04-2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check image
[im,in,ic]=size(ima);
if ic==3,
    ima = rgb2gray(ima);
end;
if isa(ima,'uint8') || isa(ima,'uint16'),
    ima=double(ima);
end;
ima=ima(:);       % vectorize ima
mi=min(ima);      % deal with negative 
ima=ima-mi+1;

% create image histogram
m=max(ima);
h=zeros(1,m);
hc=h;

for i=1:m,
    h(i)=sum(ima==i);
end;
ind=find(h);

% initiate centroids
intensityC=(1:k)*m/(k+1);

% start process
while (true),
  oldintensityC=intensityC;
  % current classification  
  for i=ind,
      [~,id]=min(abs(i-intensityC));
      hc(i)=id;
  end;
  
  %recalculation of means
  for i=1:k, 
      id=find(hc==i);
      intensityC(i)=sum(id.*h(id))/sum(h(id));
  end;
  
  if(intensityC==oldintensityC),
      break;
  end;
end;

% calculate label
% ima = reshape(ima,im,in);
% label=zeros(im,in);
% for i=1:im,
%     for j=1:in,
%         [~,id]=min(abs(ima(i,j)-intensityC));
%         label(i,j)=id;
%     end;
% end;

[~,label] = min(abs(repmat(reshape(ima,im,in),[1,1,k])-repmat(reshape(intensityC,[1,1,k]),[im,in])), [], 3);

intensityC=intensityC+mi-1;   % recover real range

