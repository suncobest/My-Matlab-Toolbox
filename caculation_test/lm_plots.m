function lm_plots ( n, t, y_dat, y_fit, sigma_y, cvg_hst )
% lm_plots ( n, t, y_dat, y_fit, sigma_y, cvg_hst )
% Plot statistics of the results of a Levenberg-Marquardt least squares
% analysis with lm.m

y_dat = y_dat(:);
y_fit = y_fit(:);

figure(100); % ------------ plot convergence history of fit
 clf
 subplot(211)
  plot( cvg_hst(:,1), cvg_hst(:,2:n+1), '-o','linewidth',4);
   ylabel('parameter values')
 
 subplot(212)
  semilogy( cvg_hst(:,1) , [ cvg_hst(:,n+2) cvg_hst(:,n+3) ], '-o','linewidth',4)
   legend('\chi^2','\lambda');
   xlabel('function calls'); ylabel('\chi^2 and \lambda')

figure(101); % ------------ plot data, fit, and confidence interval of fit
 clf
 subplot(211)
   plot(t,y_dat,'o', t,y_fit,'-', 'linewidth',2, ...
        t,y_fit+1.96*sigma_y,'-r', t,y_fit-1.96*sigma_y,'-r');
    legend('y_{data}','y_{fit}','y_{fit}+1.96\sigma_y','y_{fit}-1.96\sigma_y');
    ylabel('y(t)')
 subplot(212)
   semilogy(t,sigma_y,'-r','linewidth',4);
    xlabel('t'); ylabel('\sigma_y(t)')

figure(102); % ------------ plot histogram of residuals, are they Gaussean?
 clf
 hist(y_dat - y_fit)
  title('histogram of residuals')
  axis('tight'); xlabel('y_{data} - y_{fit}'); ylabel('count')
