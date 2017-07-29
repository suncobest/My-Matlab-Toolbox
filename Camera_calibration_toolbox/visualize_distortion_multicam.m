% visualize_distortion_multicam
%
% plot distortion models
%
% This is a slightly modified version of the script plot_CCT_distortion.m written by Mr. Oshel
% Thank you Mr. Oshel for your contribution!

if ~exist('fc_mat','var')||~exist('cc_mat','var')||~exist('kc_mat','var')||~exist('alpha_vec','var'),
    fprintf(1,'No intrinsic camera parameters available.\n');
    return;
end;

flag = 0;
if exist('fc_mat_error', 'var') && exist('err_std', 'var'),
    flag = 1;
end;

for pp = 1:n_cam,
    fc = fc_mat(:,pp);
    cc = cc_mat(:,pp);
    kc = kc_mat(:,pp);
    alpha_c = alpha_vec(pp);
    nx = imsize(1,pp);
    ny = imsize(2,pp);
    if norm(kc)==0,
        fprintf(1,'The lens distortion of camera %d seems to be zero!\n',pp);
        continue;
    end;
    fprintf(1,'Display the image distortion of camera %d:\n',pp);
    
    [gx,gy] = meshgrid(0:nx/20:(nx-1),0:ny/20:(ny-1));
    [ngx,ngy]=size(gx);
    npts = ngx*ngy;
    px=reshape(gx,1,npts);
    py=reshape(gy,1,npts);
    kk_new=[fc(1) alpha_c*fc(1) cc(1);0 fc(2) cc(2);0 0 1];
    rays=kk_new\[px;py;ones(1,length(px))];
    x=[rays(1,:)./rays(3,:);rays(2,:)./rays(3,:)];
    
    fh1 = 2;
    figure(fh1); clf;
    xd=apply_distortion(x,kc);
    px2=fc(1)*(xd(1,:)+alpha_c*xd(2,:))+cc(1);
    py2=fc(2)*xd(2,:)+cc(2);
    dx=px2-px;
    dy=py2-py;
    quiver(px+1,py+1,dx,dy);
    hold on;
    plot(cc(1)+1,cc(2)+1,'o');
    plot((nx-1)/2+1,(ny-1)/2+1,'x');
    dr=reshape(sqrt((dx.*dx)+(dy.*dy)),ngx,ngy);
    [C,h]=contour(gx,gy,dr,'k');
    clabel(C,h);
    title(['Complete Distortion Model of Camera ', num2str(pp)]);
    
    axis ij;    % "matrix" axes mode.
    axis([0 nx 0 ny]+0.5);
    axis equal;
    axis tight;
    
    if flag,
        fc_error = fc_mat_error(:,pp);
        cc_error = cc_mat_error(:,pp);
        kc_error = kc_mat_error(:,pp);
        alpha_c_error = alpha_vec_error(pp);
        
        position=get(gca,'Position');
        shr = 0.9;
        position(1)=position(1)+position(3)*((1-shr)/2);
        position(2)=position(2)+position(4)*(1-shr)+0.03;
        position(3:4)=position(3:4)*shr;
        set(gca,'position',position);
        set(gca,'fontsize',8,'fontname','clean');
        
        line1=sprintf('Principal Point               = (%0.6g, %0.6g)',cc(1),cc(2));
        line2=sprintf('Focal Length                  = (%0.6g, %0.6g)',fc(1),fc(2));
        line3=sprintf('Radial coefficients           = (%0.4g, %0.4g, %0.4g)',kc(1),kc(2),kc(5));
        line4=sprintf('Tangential coefficients       = (%0.4g, %0.4g)',kc(3),kc(4));
        line5=sprintf('+/- [%0.4g, %0.4g]',cc_error(1),cc_error(2));
        line6=sprintf('+/- [%0.4g, %0.4g]',fc_error(1),fc_error(2));
        line7=sprintf('+/- [%0.4g, %0.4g, %0.4g]',kc_error(1),kc_error(2),kc_error(5));
        line8=sprintf('+/- [%0.4g, %0.4g]',kc_error(3),kc_error(4));
        line9=sprintf('Pixel error                   = [%0.4g, %0.4g]',err_std(1),err_std(2));
        line10=sprintf('Skew                          = %0.4g',alpha_c);
        line11=sprintf('+/- %0.4g',alpha_c_error);
        
        axes('position',[0 0 1 1],'visible','off');
        text(0.11,0,{line9,line2,line1,line10,line3,line4},'horizontalalignment','left','verticalalignment','bottom','fontsize',8,'fontname','clean');
        text(0.9,0.,{line6,line5,line11,line7,line8},'horizontalalignment','right','verticalalignment','bottom','fontsize',8,'fontname','clean');
    end;
    set(fh1,'color',[1,1,1]);
    hold off;
    
    fh2 = 3;
    figure(fh2); clf;
    xd=apply_distortion(x,[0 0 kc(3) kc(4) 0]);
    px2=fc(1)*(xd(1,:)+alpha_c*xd(2,:))+cc(1);
    py2=fc(2)*xd(2,:)+cc(2);
    dx=px2-px;
    dy=py2-py;
    quiver(px+1,py+1,dx,dy);
    hold on;
    plot(cc(1)+1,cc(2)+1,'o');
    plot((nx-1)/2+1,(ny-1)/2+1,'x');
    dr=reshape(sqrt((dx.*dx)+(dy.*dy)),ngx,ngy);
    [C,h]=contour(gx,gy,dr,'k');
    clabel(C,h);
    title(['Tangential Distortion Model of Camera ', num2str(pp)]);
    
    axis ij;
    axis([0 nx 0 ny]+0.5);
    axis equal;
    axis tight;
    
    if flag,
        position=get(gca,'Position');
        shr = 0.9;
        position(1)=position(1)+position(3)*((1-shr)/2);
        position(2)=position(2)+position(4)*(1-shr)+0.03;
        position(3:4)=position(3:4)*shr;
        set(gca,'position',position);
        set(gca,'fontsize',8,'fontname','clean');
        
        line1=sprintf('Principal Point               = (%0.6g, %0.6g)',cc(1),cc(2));
        line2=sprintf('Focal Length                  = (%0.6g, %0.6g)',fc(1),fc(2));
        line3=sprintf('Radial coefficients           = (%0.4g, %0.4g, %0.4g)',kc(1),kc(2),kc(5));
        line4=sprintf('Tangential coefficients       = (%0.4g, %0.4g)',kc(3),kc(4));
        line5=sprintf('+/- [%0.4g, %0.4g]',cc_error(1),cc_error(2));
        line6=sprintf('+/- [%0.4g, %0.4g]',fc_error(1),fc_error(2));
        line7=sprintf('+/- [%0.4g, %0.4g, %0.4g]',kc_error(1),kc_error(2),kc_error(5));
        line8=sprintf('+/- [%0.4g, %0.4g]',kc_error(3),kc_error(4));
        line9=sprintf('Pixel error                   = [%0.4g, %0.4g]',err_std(1),err_std(2));
        line10=sprintf('Skew                          = %0.4g',alpha_c);
        line11=sprintf('+/- %0.4g',alpha_c_error);
        
        axes('position',[0 0 1 1],'visible','off');
        text(0.11,0,{line9,line2,line1,line10,line3,line4},'horizontalalignment','left','verticalalignment','bottom','fontsize',8,'fontname','clean');
        text(0.9,0.,{line6,line5,line11,line7,line8},'horizontalalignment','right','verticalalignment','bottom','fontsize',8,'fontname','clean');
    end;
    set(fh2,'color',[1,1,1]);
    hold off;
    
    
    fh3 = 4;
    figure(fh3); clf;
    xd=apply_distortion(x,[kc(1) kc(2) 0 0 kc(5)]);
    px2=fc(1)*(xd(1,:)+alpha_c*xd(2,:))+cc(1);
    py2=fc(2)*xd(2,:)+cc(2);
    dx=px2-px;
    dy=py2-py;
    quiver(px+1,py+1,dx,dy);
    hold on;
    plot(cc(1)+1,cc(2)+1,'o');
    plot((nx-1)/2+1,(ny-1)/2+1,'x');
    dr=reshape(sqrt((dx.*dx)+(dy.*dy)),ngx,ngy);
    [C,h]=contour(gx,gy,dr,'k');
    clabel(C,h);
    title(['Radial Distortion Model of Camera ', num2str(pp)]);
    
    axis ij;
    axis([0 nx 0 ny]+0.5);
    axis equal;
    axis tight;
    
    if flag,
        position=get(gca,'Position');
        shr = 0.9;
        position(1)=position(1)+position(3)*((1-shr)/2);
        position(2)=position(2)+position(4)*(1-shr)+0.03;
        position(3:4)=position(3:4)*shr;
        set(gca,'position',position);
        set(gca,'fontsize',8,'fontname','clean');
        
        line1=sprintf('Principal Point               = (%0.6g, %0.6g)',cc(1),cc(2));
        line2=sprintf('Focal Length                  = (%0.6g, %0.6g)',fc(1),fc(2));
        line3=sprintf('Radial coefficients           = (%0.4g, %0.4g, %0.4g)',kc(1),kc(2),kc(5));
        line4=sprintf('Tangential coefficients       = (%0.4g, %0.4g)',kc(3),kc(4));
        line5=sprintf('+/- [%0.4g, %0.4g]',cc_error(1),cc_error(2));
        line6=sprintf('+/- [%0.4g, %0.4g]',fc_error(1),fc_error(2));
        line7=sprintf('+/- [%0.4g, %0.4g, %0.4g]',kc_error(1),kc_error(2),kc_error(5));
        line8=sprintf('+/- [%0.4g, %0.4g]',kc_error(3),kc_error(4));
        line9=sprintf('Pixel error                   = [%0.4g, %0.4g]',err_std(1),err_std(2));
        line10=sprintf('Skew                          = %0.4g',alpha_c);
        line11=sprintf('+/- %0.4g',alpha_c_error);
        
        axes('position',[0 0 1 1],'visible','off');
        text(0.11,0,{line9,line2,line1,line10,line3,line4},'horizontalalignment','left','verticalalignment','bottom','fontsize',8,'fontname','clean');
        text(0.9,0.,{line6,line5,line11,line7,line8},'horizontalalignment','right','verticalalignment','bottom','fontsize',8,'fontname','clean');
    end;
    set(fh3,'color',[1,1,1]);
    hold off;
    figure(fh1);
    
    fprintf(1,'\nPlease check the distortion models of camera %d.\n', pp);
    if pp<n_cam,
        fprintf(1, 'press any key to continue ...\n');
        pause;
    end;
end;

disp('Done');