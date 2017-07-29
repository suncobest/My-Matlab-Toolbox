% 列出12种欧拉角惯例，如XZX，YZX...

% 统一表示绕各轴的基本旋转矩阵
%     Rx=[1,0,0; 0,cos(ax),-sin(ax); 0,sin(ax),cos(ax)];
%     Ry=[cos(ay),0,sin(ay); 0,1,0; -sin(ay),0,cos(ay)];
%     Rz=[cos(az),-sin(az),0; sin(az),cos(az),0; 0,0,1];

% pwd返回当前工作路径，genpath(dir)返回所有dir文件夹以内的所有子文件夹路径。
addpath(genpath(pwd));  
clear; clc; echo off;
syms a1 a2 a3  % 欧拉角变量，绕3个轴的转角依次为：a1，a2，a3

name='XYZ';
angle = ['ax';'ay';'az'];
% R = {'[1,0,0; 0,cos(ax),-sin(ax); 0,sin(ax),cos(ax)]',...
%     '[cos(ay),0,sin(ay); 0,1,0; -sin(ay),0,cos(ay)]',...
%     '[cos(az),-sin(az),0; sin(az),cos(az),0; 0,0,1]'};

Rx = '[1,0,0; 0,cos(ax),-sin(ax); 0,sin(ax),cos(ax)]';
Ry = '[cos(ay),0,sin(ay); 0,1,0; -sin(ay),0,cos(ay)]';
Rz = '[cos(az),-sin(az),0; sin(az),cos(az),0; 0,0,1]';
% char(Rx,Ry,Rz)将文本字符串Rx，Ry，Rz组成字符矩阵，每个字符串作为一行；
% 以字符最多的行作为基准，其他行在末尾补充空格，以保证矩阵成立。
R = char(Rx,Ry,Rz); 

idx = perm(1:3,2);
idx = [idx, idx(:,1)];
idx(7:12,:) = perm(1:3);  % idx的每行是旋转轴的次序，如[1 2 1],[1 2 3],...

% XYZ = perm('XYZ',2);  % 从'XYZ'中选出2个字母取全排列（6行2列），每个字符占用一列。
% XYZ = [XYZ, XYZ(:,1)];  % 前6行是6种严格意义上的欧拉角，有两根相同的转轴，如XYX，ZXZ...
% XYZ(7:12,:) = perm('XYZ');  % 后6行是6种Tait-Bryan角，如XYZ，YZX...

fid = fopen('Euler_angles_conventions.txt','w+');
% fid = fopen('Euler_angles_conventions.m','w+');
fprintf(fid,'%%%% 列出12种欧拉角惯例\r\n\r\n');
fprintf(fid,'clear all\r\n');
fprintf(fid,'syms a1 a2 a3\r\n\r\n');

for i=1:12
    A=angle(idx(i,:),:);  % 利用行指标对angle进行“行操作”，置换或复制
    % 由于欧拉角变量名为ax，ay，az；对于严格意义的欧拉角，如ZXZ，第一次和第三次都是绕Z轴旋转，需要对az进行两次赋值；
    % 所以首先对前两个欧拉角赋值，az=a1，ax=a2，计算前两次的旋转矩阵ZX；然后再对第三个欧拉角赋值，az=a3，计算ZXZ。
    eval([A(1,:),'=a1']);   
    eval([A(2,:),'=a2']);   % 给前两个欧拉角赋值；A(1,:)和A(2,:)为字符串变量
    N=name(idx(i,:));
%     rotation=eval([R{idx(i,1)},'*',R{idx(i,2)}]);  % R{idx(i,1)}为元胞中的字符串，eval(R{idx(i,1)})的含义是矩阵
    
    rotation = eval([R(idx(i,1),:), '*', R(idx(i,2),:)]);
    
    eval([A(3,:),'=a3']);   % 给第3个欧拉角赋值
%     rotation=eval(['rotation*', R{idx(i,3)}]);  % 根据XYZ=Rx*Ry*Rz, ZXZ=Rz*Rx*Rz,...
    
    rotation = rotation * eval( R(idx(i,3),:));
    eval([N,'=rotation']);
    
    fprintf(fid,'%%%% --------------------------------------------------------------------------------------------\r\n\r\n');


%    下面两种书写方式等价，无所谓txt或m文件。    
%%  写txt文件 
    for a=1:3
        fprintf(fid,'%3s(%d, :) = [ ', N, a);
        for b=1:2
            s=char(rotation(a,b));  % rotation是对于符号变量
            fprintf(fid,'%s,   ',s);
        end
        s=char(rotation(a,3));
        fprintf(fid,'%s ];\r\n',s);
    end
   
    
%%  写m文件    
%     str = char(rotation);
%     % 由于rotation是sym符号变量，上式会返回类似matrix([[rotation(1,:)], [rotation(2,:)], [rotation(3,:)]])
%     % 形式的一行字符串。
%     str(1:length('matrix([')) = [];  % 删除首尾多余的字符
%     str(end-1:end) = [];
%     id1 = find(str=='[');   % 找出rotation每一行的开始位置的指标
%     id2 = find(str==']');   % 找出rotation每一行的结束位置的指标
%     for j=1:3
%         fprintf( fid,'%3s(%d, :) = %s;\r\n', N, j, str(id1(j):id2(j)) );
%     end


%%    
    
    fprintf(fid,'\r\n\r\n');
    
end

fclose(fid);
