% Curve fitting with the steepest gradient descent algorithm
% Problem: Fit a set of point data {xi,yi} (xi=0, 0.1, ..., 2*PI) with the function: 
% f(x) = a*cos(bx)+b*sin(ax)
% where a and b are the parameters that we want to estimate by the LM method.

% 1. Using the MATLAB symbolic toolbox, find the analytic form of the Jacobian of one
% row of d with respect to the parameters a and b

% syms a b x y real;
% f=( a * cos(b*x) + b * sin(a*x));
% d=y-f;
% Jsym=jacobian(d,[a;b]);
% 
% % This MATLAB code returns the following symbolic expression
% Jsym =[ -cos(b*x)-b*cos(a*x)*x, a*sin(b*x)*x-sin(a*x)]
% Jp=jacobian(Jsym,[a,b])
% Jp=[b*x^2*sin(a*x), x*sin(b*x) - x*cos(a*x); x*sin(b*x) - x*cos(a*x), a*x^2*cos(b*x)];
% 2. Generate the synthetic data from the curve function with some additional noise

a=100;
b=102;
x=(0:0.1:2*pi)';
nx=length(x);
y = a * cos(b*x) + b * sin(a*x);
% add random noise
y_input = y + 8*randn(nx,1);


% 3.The main code for the SG algorithm for estimating a and b from the above data

% initial guess for the parameters
a0=100.5; b0=101.4;
y_init = a0 * cos(b0*x) + b0 * sin(a0*x);
Ndata=length(y_input);
Nparams=2; % a and b are the parameters to be estimated
n_iters=100; % set # of iterations for the LM
a_est=a0;
b_est=b0;

for it=1:n_iters
    y_est = a_est * cos(b_est*x) + b_est * sin(a_est*x);
    d=y_input-y_est;
    e=dot(d,d)
    % Evaluate the Jacobian matrix at the current parameters (a_est, b_est)
    J=[-cos(b_est*x)-(b_est*cos(a_est*x).*x), (a_est*sin(b_est*x).*x)-sin(a_est*x)];
    gradient=J'*d(:);
    % compute the approximated Hessian matrix, J is the transpose of J
    H=J'*J;
    % line search, step length for steepest gradient descent(zig-zag)
    lamda=(gradient'*gradient)/(gradient'*H*gradient);
    dp=-lamda*gradient;
    a_est=a_est+dp(1);
    b_est=b_est+dp(2);
    
    if norm(gradient)<1e-7
        break;
    end
    createfigure(x, [y_input,y_est]);
    pause(0.2)
end




