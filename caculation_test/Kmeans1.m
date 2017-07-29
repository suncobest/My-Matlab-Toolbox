a=5*[randn(20,1)-5,randn(20,1)+5];
b=5*[randn(20,1)+5,randn(20,1)+5];
c=5*[randn(20,1),randn(20,1)];
d=5*[randn(20,1)-5,randn(20,1)-5];
e=5*[randn(20,1)+5,randn(20,1)-5];
Data=[a;b;c;d;e];
k=5;
[IDX,C,SUMD,D]=kmeans(Data,k);
% % plot(a(:,1),a(:,2),'.');hold on;
% % plot(b(:,1),b(:,2),'r.');hold on;
% % plot(c(:,1),c(:,2),'c.');hold on;
% % plot(d(:,1),d(:,2),'g.');hold on;
% % plot(e(:,1),e(:,2),'y.');hold on;

% I1=find(IDX==1);
% I2=find(IDX==2);
% I3=find(IDX==3);
% I4=find(IDX==4);
% I5=find(IDX==5);
% plot(Data(I1,1),Data(I1,2),'b.',C(1,1),C(1,2),'b*','MarkerSize',9);hold on;
% plot(Data(I2,1),Data(I2,2),'r.',C(2,1),C(2,2),'r*','MarkerSize',9);hold on;
% plot(Data(I3,1),Data(I3,2),'g.',C(3,1),C(3,2),'g*','MarkerSize',9);hold on;
% plot(Data(I4,1),Data(I4,2),'m.',C(4,1),C(4,2),'m*','MarkerSize',9);hold on;
% plot(Data(I5,1),Data(I5,2),'y.',C(5,1),C(5,2),'y*','MarkerSize',9);hold on;
% for i=1:100
%     text(Data(i,1),Data(i,2),num2str(IDX(i)))
% end
color='brgmy';
I=cell(1,k);
figure(1);hold on;
for j=1:k
    I{1,j}=find(IDX==j);
    plot(Data(I{1,j},1),Data(I{1,j},2),[color(j),'.'],...
        C(j,1),C(j,2),[color(j),'*'],'MarkerSize',9);
    text(Data(I{1,j},1),Data(I{1,j},2),num2str(j),'color',color(j));
end


print(gcf,'-dpng','kmeans1.png');