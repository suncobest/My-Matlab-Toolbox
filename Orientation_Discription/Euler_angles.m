% �г�12��ŷ���ǹ�������XZX��YZX...

% ͳһ��ʾ�Ƹ���Ļ�����ת����
%     Rx=[1,0,0; 0,cos(ax),-sin(ax); 0,sin(ax),cos(ax)];
%     Ry=[cos(ay),0,sin(ay); 0,1,0; -sin(ay),0,cos(ay)];
%     Rz=[cos(az),-sin(az),0; sin(az),cos(az),0; 0,0,1];

% pwd���ص�ǰ����·����genpath(dir)��������dir�ļ������ڵ��������ļ���·����
addpath(genpath(pwd));  
clear; clc; echo off;
syms a1 a2 a3  % ŷ���Ǳ�������3�����ת������Ϊ��a1��a2��a3

name='XYZ';
angle = ['ax';'ay';'az'];
% R = {'[1,0,0; 0,cos(ax),-sin(ax); 0,sin(ax),cos(ax)]',...
%     '[cos(ay),0,sin(ay); 0,1,0; -sin(ay),0,cos(ay)]',...
%     '[cos(az),-sin(az),0; sin(az),cos(az),0; 0,0,1]'};

Rx = '[1,0,0; 0,cos(ax),-sin(ax); 0,sin(ax),cos(ax)]';
Ry = '[cos(ay),0,sin(ay); 0,1,0; -sin(ay),0,cos(ay)]';
Rz = '[cos(az),-sin(az),0; sin(az),cos(az),0; 0,0,1]';
% char(Rx,Ry,Rz)���ı��ַ���Rx��Ry��Rz����ַ�����ÿ���ַ�����Ϊһ�У�
% ���ַ���������Ϊ��׼����������ĩβ����ո��Ա�֤���������
R = char(Rx,Ry,Rz); 

idx = perm(1:3,2);
idx = [idx, idx(:,1)];
idx(7:12,:) = perm(1:3);  % idx��ÿ������ת��Ĵ�����[1 2 1],[1 2 3],...

% XYZ = perm('XYZ',2);  % ��'XYZ'��ѡ��2����ĸȡȫ���У�6��2�У���ÿ���ַ�ռ��һ�С�
% XYZ = [XYZ, XYZ(:,1)];  % ǰ6����6���ϸ������ϵ�ŷ���ǣ���������ͬ��ת�ᣬ��XYX��ZXZ...
% XYZ(7:12,:) = perm('XYZ');  % ��6����6��Tait-Bryan�ǣ���XYZ��YZX...

fid = fopen('Euler_angles_conventions.txt','w+');
% fid = fopen('Euler_angles_conventions.m','w+');
fprintf(fid,'%%%% �г�12��ŷ���ǹ���\r\n\r\n');
fprintf(fid,'clear all\r\n');
fprintf(fid,'syms a1 a2 a3\r\n\r\n');

for i=1:12
    A=angle(idx(i,:),:);  % ������ָ���angle���С��в��������û�����
    % ����ŷ���Ǳ�����Ϊax��ay��az�������ϸ������ŷ���ǣ���ZXZ����һ�κ͵����ζ�����Z����ת����Ҫ��az�������θ�ֵ��
    % �������ȶ�ǰ����ŷ���Ǹ�ֵ��az=a1��ax=a2������ǰ���ε���ת����ZX��Ȼ���ٶԵ�����ŷ���Ǹ�ֵ��az=a3������ZXZ��
    eval([A(1,:),'=a1']);   
    eval([A(2,:),'=a2']);   % ��ǰ����ŷ���Ǹ�ֵ��A(1,:)��A(2,:)Ϊ�ַ�������
    N=name(idx(i,:));
%     rotation=eval([R{idx(i,1)},'*',R{idx(i,2)}]);  % R{idx(i,1)}ΪԪ���е��ַ�����eval(R{idx(i,1)})�ĺ����Ǿ���
    
    rotation = eval([R(idx(i,1),:), '*', R(idx(i,2),:)]);
    
    eval([A(3,:),'=a3']);   % ����3��ŷ���Ǹ�ֵ
%     rotation=eval(['rotation*', R{idx(i,3)}]);  % ����XYZ=Rx*Ry*Rz, ZXZ=Rz*Rx*Rz,...
    
    rotation = rotation * eval( R(idx(i,3),:));
    eval([N,'=rotation']);
    
    fprintf(fid,'%%%% --------------------------------------------------------------------------------------------\r\n\r\n');


%    ����������д��ʽ�ȼۣ�����νtxt��m�ļ���    
%%  дtxt�ļ� 
    for a=1:3
        fprintf(fid,'%3s(%d, :) = [ ', N, a);
        for b=1:2
            s=char(rotation(a,b));  % rotation�Ƕ��ڷ��ű���
            fprintf(fid,'%s,   ',s);
        end
        s=char(rotation(a,3));
        fprintf(fid,'%s ];\r\n',s);
    end
   
    
%%  дm�ļ�    
%     str = char(rotation);
%     % ����rotation��sym���ű�������ʽ�᷵������matrix([[rotation(1,:)], [rotation(2,:)], [rotation(3,:)]])
%     % ��ʽ��һ���ַ�����
%     str(1:length('matrix([')) = [];  % ɾ����β������ַ�
%     str(end-1:end) = [];
%     id1 = find(str=='[');   % �ҳ�rotationÿһ�еĿ�ʼλ�õ�ָ��
%     id2 = find(str==']');   % �ҳ�rotationÿһ�еĽ���λ�õ�ָ��
%     for j=1:3
%         fprintf( fid,'%3s(%d, :) = %s;\r\n', N, j, str(id1(j):id2(j)) );
%     end


%%    
    
    fprintf(fid,'\r\n\r\n');
    
end

fclose(fid);
