function [XL,XR] = stereo_triangulation(xL,xR,om,T,fcl,ccl,kcl,alpha_l,fcr,ccr,kcr,alpha_r)

% [XL,XR] = stereo_triangulation(xL,xR,om,T,fcl,ccl,kcl,alpha_l,fcr,ccr,kcr,alpha_r),
%
% Function that computes the position of a set on N points given the left and right image projections.
% The cameras are assumed to be calibrated, intrinsically, and extrinsically.
%
% Input:
%           xL: 2xN matrix of pixel coordinates in the left image
%           xR: 2xN matrix of pixel coordinates in the right image
%           om,T: rotation vector and translation vector between right and left cameras (output of stereo calibration)
%           fcl,ccl,...: intrinsic parameters of the left camera  (output of stereo calibration)
%           fcr,ccr,...: intrinsic parameters of the right camera (output of stereo calibration)
%
% Output:
%
%           XL: 3xN matrix of coordinates of the points in the left camera reference frame
%           XR: 3xN matrix of coordinates of the points in the right camera reference frame
%
% Note: XR and XL are related to each other through the rigid motion equation: XR = R * XL + T, where R = rodrigues(om)
% For more information, visit http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example5.html
%
%
% (c) Jean-Yves Bouguet - Intel Corporation - April 9th, 2003



%--- Normalize the image projection according to the intrinsic parameters of the left and right cameras
xnl = normalize_pixel(xL,fcl,ccl,kcl,alpha_l);
xnr = normalize_pixel(xR,fcr,ccr,kcr,alpha_r);

%--- Number of points:
N = size(xnl,2);

%--- Extend the normalized projections in homogeneous coordinates
xnl = [xnl;ones(1,N)];
xnr = [xnr;ones(1,N)];

%--- Rotation matrix corresponding to the rigid motion between left and right cameras:
R = rodrigues(om);


%--- Triangulation of the rays in 3D space:
%  XR = R * XL + T������XR��R * XL��T�������ʸ�������ε����ߣ�������ϵCR��ͶӰ��
%  ����xnl��xnr�ֱ�ΪXL��XR�Ĺ�һ�����꣨xnl=XL/Zl,xnr=XR/Zr��������xnl��xnr�ֱ�ƽ����XL��XR��
%  ��֪xnl��xnr����֪����XR��R * XL�ķ��򣬻���֪��һ��T�����Ǽ�һ�ߣ�ʸ�������ζ��⡣


u = R * xnl;  % ��ʸ��xnl������ת��������ϵCL�е�ʸ��ת��Ϊ����ϵCR�е�ʸ����

% ��a=norm(xnl),b=norm(xnr),c=norm(T)��a=norm(u)����a��b��c�����ߵ���Ӧ�Խ�
% �ֱ�ΪA��B��C������A+B+C=pi�������Ҷ����a*Zl/sin(A)=b*Zr/sin(B)=c/sin(C)��
% ����Zl=c*sin(A)/(a*sin(C));Zr=c*sin(B)/(b*sin(C))��

xnl2 = dot(xnl,xnl);   % xnl2 = a^2,
xnr2 = dot(xnr,xnr);  % xnr2 = b^2

T_vect = repmat(T, [1 N]);

DD = xnl2 .* xnr2 - dot(u,xnr).^2;  % DD=a^2*b^2*sin(C)^2

dot_uT = dot(u,T_vect);  % dot_uT = -a*c*cos(B)
dot_xnrT = dot(xnr,T_vect);    % dot_xnrT = b*c*cos(A)
dot_xnru = dot(u,xnr);  % dot_xnru = a*b*cos(C)

NN1 = dot_xnru.*dot_xnrT - xnr2 .* dot_uT;  % NN1 = a*b^2*c*sin(A)*sin(C)
NN2 = xnl2.*dot_xnrT - dot_uT.*dot_xnru;  % NN2 = a^2*b*c*sin(B)*sin(C)

Zl = NN1./DD;
Zr = NN2./DD;

X1 = xnl .* repmat(Zl,[3 1]);
X2 = R'*(xnr.*repmat(Zr,[3,1])  - T_vect);  % ԭ����X1=X2�����Ƕ�Ӧ�㣬����ȣ�


%--- Left coordinates:
XL = 1/2 * (X1 + X2);   % ȡƽ��ֵ��X1��X2���е㣩

%--- Right coordinates:
XR = R*XL + T_vect;

