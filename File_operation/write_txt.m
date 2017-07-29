%%%  Ð´mÎÄ¼þ 
%% homography_jacobi
fid=fopen('homography_jacobi_test.m','a+');  % 'r+','w+','a+'
[m,n] = size(J);

for i=1:m
    for j=1:n
        fprintf(fid,'J(%d, %d) = %s;\r\n',i,j,char(J(i,j)) );
    end
end

fclose(fid);


%% dRdom
fid=fopen('rodrigues_jacobi.m','a+');  % 'r+','w+','a+'
[m,n] = size(dRdom);

for i=1:m
    for j=1:n
        fprintf(fid,'dRdom(%d, %d) = %s;\r\n',i,j,char(dRdom(i,j)) );
    end
end

fclose(fid);

%% domdR
fid=fopen('rodrigues_jacobi.m','a+');  % 'r+','w+','a+'
[m,n] = size(domdR);

for i=1:m
    for j=1:n
        fprintf(fid,'domdR(%d, %d) = %s;\r\n',i,j,char(domdR(i,j)) );
    end
end

fclose(fid);


%%  JA JB
fid=fopen('cam_jacobi_nodist_test.m','a+');   %('cam_jacobi_test.m','a+');  % 'r+','w+','a+'
[ma,na] = size(A);
[mb,nb] = size(B);

% for i=1:2
%     fprintf(fid,'xd(%d,1) = %s;\r\n',i, char(xd(i)) );
% end

fprintf(fid,'\r\n');
for i=1:2
    fprintf(fid,'xn(%d,1) = %s;\r\n',i, char(xn(i)) );
end

fprintf(fid,'\r\n');
for i=1:3
    fprintf(fid,'Xc(%d,1) = %s;\r\n',i, char(Xc(i)) );
end

fprintf(fid,'\r\n');
for i=1:3
    for j=1:3
        fprintf(fid,'Rc(%d, %d) = %s;\r\n',i,j,char(Rc(i,j)) );
    end
end

fprintf(fid,'\r\n\r\n');

for i=1:ma
    for j=1:na
        fprintf(fid,'JA(%d, %d) = %s;\r\n',i,j,char(A(i,j)) );
    end
end

fprintf(fid,'\r\n\r\n');

for i=1:mb
    for j=1:nb
        fprintf(fid,'JB(%d, %d) = %s;\r\n',i,j,char(B(i,j)) );
    end
end



fclose(fid);

%% write txt

fid=fopen('border2.txt','w+');  % 'r+','w+','a+'
% fprintf(fid,'x      y      z\r\n');
a = border';  % border is size of (3,N)
[m n]=size(a);
for i = 1:m
    fprintf(fid,'%12.8f        %12.8f        %12.8f\r\n',a(i,1),a(i,2),a(i,3));
end
fclose(fid);






%%
fid=fopen('Smoothed5.txt','wt+');  % 'r+','w+','a+'
fprintf(fid,'%10.8f     %10.8f\n',x1');
fclose(fid);