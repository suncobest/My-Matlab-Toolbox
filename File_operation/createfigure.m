function createfigure(X1, YMatrix1)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 22-Jun-2014 10:50:29

% Create figure
set(gcf,'Name','Curve fitting test1','Color',[1 1 1]);

% Create axes
axes1 = axes('Parent',gcf);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1);
set(plot1(1),'LineStyle',':','DisplayName','ymeas');
set(plot1(2),'Color',[1 0 0],'DisplayName','yest');

axis([0 6.5 -250 250]);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.15 0.81 0.15 0.1]);



