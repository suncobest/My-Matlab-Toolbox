function ns = count_squares(I,x1,y1,x2,y2,win)

[ny,nx] = size(I);

if ((x1-win <= 0) || (x1+win >= nx) || (y1-win <= 0) || (y1+win >= ny) || ...
        (x2-win <= 0) || (x2+win >= nx) || (y2-win <= 0) || (y2+win >= ny))
    ns = -1;
    return;
end;

if ((x1 - x2)^2+(y1-y2)^2) <  win,
    ns = -1;
    return;
end;

% �˺�����Ҫ��lambda���಻�Գ�,��ʹ��out_ff�������棨����Ҫ��������������
% lambda is the line pass through homogeneous points [x1;y1;1] and [x2;y2;1]
lambda = [y1 - y2;x2 - x1;x1*y2 - x2*y1];                 % lambda = cross([x1;y1;1],[x2;y2;1])
lambda = 1/sqrt(lambda(1)^2 + lambda(2)^2) * lambda;      % ��ʽֱ�߷���[a b c]*x/sqrt(a^2+b^2)

% l1 = lambda + [0;0;win];
% l2 = lambda - [0;0;win]; 
% ֱ��l1��l2ƽ����ֱ��lambda���ֱ��������࣬���붼��win����Ȥ����Ϊl1��l2֮��Ĵ�״����

dx = x2 - x1;
dy = y2 - y1;

if abs(dx) > abs(dy),                                     % ֱ��lambda��x��н�С��45�ȣ�ȡy=f(x)������ȡx=f(y)
   if x2 > x1,
      xs = x1:x2;
   else
      xs = x1:-1:x2;
   end;
   ys = -(lambda(3) + lambda(1)*xs)/lambda(2);             % ����[xs;ys;1]��ֱ��lambda��
else
   if y2 > y1,
       ys = y1:y2;
   else
       ys = y1:-1:y2;
   end;
   xs = -(lambda(3) + lambda(2)*ys)/lambda(1);
end;

Np = length(xs);
xs_mat = ones(2*win + 1,1)*xs;
ys_mat = ones(2*win + 1,1)*ys;
win_mat = (-win:win)'*ones(1,Np);
xs_mat2 = round(xs_mat - win_mat * lambda(1));              % �������[xs,ys]-n*[a,b]��n=-win:win������norm([a,b])=1
ys_mat2 = round(ys_mat - win_mat * lambda(2));              % �õ���ֱ��lambda����-win��win�����ص�����
ind_mat = (xs_mat2 - 1) * ny + ys_mat2;                     % ����������[x,y]ת��Ϊ�����±�(x-1)*ny+y
ima_patch = zeros(2*win + 1,Np);
ima_patch(:) = I(ind_mat(:));                               % ��״��������ػҶ�

%ima2 = ima_patch(:,win+1:end-win);

filtk = [ones(win,Np);zeros(1,Np);-ones(win,Np)];           % ���ɴ�״ϵ�����ϰ벿�֣���win��Ϊ1���м�һ��Ϊ0���°벿��Ϊ-1
out_f = sum(filtk.*ima_patch);                              % ����״��������ػҶ��ϰ벿�ּ�ȥ�°벿�֣�����ͣ��ش�ֱ��lambda����
out_ff = conv2(out_f,[1/4 1/2 1/4],'same');                % out_fΪ��������һά�������ʹ��conv��conv2û��������������tentģ��
out_ff = out_ff(win+1:end-win);                           % ȥ�������˵㣺����ͷ�ֱ��ȥ�뾶win�ĳ���
ns = length(find(((out_ff(2:end)>=0)&(out_ff(1:end-1)<0)) | ((out_ff(2:end)<=0)&(out_ff(1:end-1)>0))))+1; 

% Ѱ��out_ff�ı���ŵĵ�ĸ������м�Ľǵ�������Ȼ��+1�õ���������
% (out_ff(2:end)>=0)&(out_ff(1:end-1)<0)�ҳ�n-1Ϊ����nΪ�Ǹ���ת�۵㣨������-��+��
% (out_ff(2:end)<=0)&(out_ff(1:end-1)>0)�ҳ�n-1Ϊ����nΪ������ת�۵㣨������+��-��
% �˺��������ҳ���0��-��ת�۵㣬Ҳ�����ҳ���0��+��ת�۵�

return;
