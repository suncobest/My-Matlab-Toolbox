% Curve fitting with the LM algorithm
% Problem: Fit a set of point data {xi,yi} (xi=0, 0.1, ..., 2*PI) with the function: 
% f(x) = a*cos(bx)+b*sin(ax)
% where a and b are the parameters that we want to estimate by the LM method.

% 1. Using the MATLAB symbolic toolbox, find the analytic form of the Jacobian of one
% row of d with respect to the parameters a and b

% syms a b x y real;
% f=( a * cos(b*x) + b * sin(a*x));
% d=y-f;
% Jsym=jacobian(d,[a b]);

% This MATLAB code returns the following symbolic expression
% Jsym =[ -cos(b*x)-b*cos(a*x)*x, a*sin(b*x)*x-sin(a*x)]
% 
% 2. Generate the synthetic data from the curve function with some additional noise

a=100;
b=102;
x=(0:0.1:2*pi)';
y = a * cos(b*x) + b * sin(a*x);
% add random noise
y_input = y + 8*randn(length(x),1);


% 3.The main code for the LM algorithm for estimating a and b from the above data

% initial guess for the parameters
a0=99.4; b0=101.4;
y_init = a0 * cos(b0*x) + b0 * sin(a0*x);
Ndata=length(y_input);
Nparams=2; % a and b are the parameters to be estimated
n_iters=100; % set # of iterations for the LM
lamda=0.01; % set an initial value of the damping factor for the LM
updateJ=1;
a_est=a0;
b_est=b0;
y_est = a_est * cos(b_est*x) + b_est * sin(a_est*x);
d=y_input-y_est;
e=dot(d,d);  % the inital error      
for it=1:n_iters
    if updateJ==1
        % Evaluate the Jacobian matrix at the current parameters (a_est, b_est)
        J=[-cos(b_est*x)-(b_est*cos(a_est*x).*x), (a_est*sin(b_est*x).*x)-sin(a_est*x)];
        % Evaluate the distance error at the current parameters
        y_est = a_est * cos(b_est*x) + b_est * sin(a_est*x);
        d=y_input-y_est;
        % compute the approximated Hessian matrix, J�� is the transpose of J
        H=J'*J;
    end
    % Apply the damping factor to the Hessian matrix
    H_lm=H+(lamda*eye(Nparams));
    % Compute the updated parameters
    dp=-H_lm\(J'*d(:));     % dp=-inv(H_lm)*(J'*d(:));
    a_lm=a_est+dp(1);
    b_lm=b_est+dp(2);
    % Evaluate the total distance error at the updated parameters
    y_est_lm = a_lm * cos(b_lm*x) + b_lm * sin(a_lm*x);
    d_lm=y_input-y_est_lm;
    e_lm=dot(d_lm,d_lm)
    % If the total distance error of the updated parameters is less than the previous one
    % then makes the updated parameters to be the current parameters
    % and decreases the value of the damping factor
    if e_lm<e
        lamda=lamda/10;
        a_est=a_lm;
        b_est=b_lm;
        e=e_lm;
        updateJ=1;
    else % otherwise increases the value of the damping factor
        updateJ=0;
        lamda=lamda*10;
    end
    gradient = J'*d_lm(:);
    if norm(gradient)<1e-7
        break;
    end
    createfigure(x, [y_input,y_est_lm]);
    pause(0.2)
end




