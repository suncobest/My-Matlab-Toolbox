function diff_sig = anisodiff1D(sig, kappa, num_iter, flag, delta_t)
%ANISODIFF1D Conventional anisotropic diffusion
%   DIFF_SIG = ANISODIFF1D(SIG, KAPPA, NUM_ITER, flag, DELTA_T) perfoms 
%   conventional anisotropic diffusion (Perona & Malik) upon 1D signals.
%   A 1D network structure of 2 neighboring nodes is considered for diffusion 
%   conduction.
% 
%       ARGUMENT DESCRIPTION:
%               SIG      - input column vector (Nx1).
%               KAPPA    - gradient modulus threshold that controls the conduction.
%               NUM_ITER - number of iterations. 
%               flag   - switch of conduction coefficient functions proposed by Perona & Malik:
%                          0 - c(x,t) = exp(-(nablaI/kappa).^2),
%                              privileges high-contrast edges over low-contrast ones. 
%                          1 - c(x,t) = 1./(1 + (nablaI/kappa).^2),
%                              privileges wide regions over smaller ones.  (Default)
%               DELTA_T  - integration constant (0 <= delta_t <= 1/3).
%                          Usually, due to numerical stability this 
%                          parameter is set to its maximum value.
%
%       OUTPUT DESCRIPTION:
%               DIFF_SIG - (diffused) signal with the largest scale-space parameter.
% 
%   Example
%   -------------
%   s = exp(-(-10:0.1:10).^2/2.5)' + rand(length(-10:0.1:10),1);
%   num_iter = 15;
%   kappa = 30;
%   ad = anisodiff1D(s, kappa, num_iter);
%   figure, subplot 121, plot(s), axis([0 200 0 2]), subplot 122, plot(ad), axis([0 200 0 2])
% 
% See also anisodiff2D, anisodiff3D.

% References: 
%   P. Perona and J. Malik. 
%   Scale-Space and Edge Detection Using Anisotropic Diffusion.
%   IEEE Transactions on Pattern Analysis and Machine Intelligence, 
%   12(7):629-639, July 1990.
% 
%   G. Grieg, O. Kubler, R. Kikinis, and F. A. Jolesz.
%   Nonlinear Anisotropic Filtering of MRI Data.
%   IEEE Transactions on Medical Imaging,
%   11(2):221-232, June 1992.
% 
%   MATLAB implementation based on Peter Kovesi's anisodiff(.):
%   P. D. Kovesi. MATLAB and Octave Functions for Computer Vision and Image Processing.
%   School of Computer Science & Software Engineering,
%   The University of Western Australia. Available from:
%   <http://www.csse.uwa.edu.au/~pk/research/matlabfns/>.
% 
% Credits:
% Daniel Simoes Lopes
% ICIST
% Instituto Superior Tecnico - Universidade Tecnica de Lisboa
% danlopes (at) civil ist utl pt
% http://www.civil.ist.utl.pt/~danlopes
%
% May 2007 original version.

if nargin<5,
    delta_t = 1/3;
    if nargin<4,
        flag = 1;
        if nargin<3,
            num_iter = 15;
            if nargin<2,
                kappa = 10;
            end;
        end;
    end;
end;
flag = ~~flag;

% Convert input signal to double.
sig = double(sig);

% PDE (partial differential equation) initial condition.
diff_sig = sig;

% Center point distance.
dx = 1;

% 1D convolution masks - finite differences.
hW = [1 -1 0]';
hE = [0 -1 1]';

% Anisotropic diffusion.
for t = 1:num_iter,

        % Finite differences. [imfilter(.,.,'conv') can be replaced by conv(.,.,'same')]
        nablaW = imfilter(diff_sig,hW,'conv');
        nablaE = imfilter(diff_sig,hE,'conv');   

        % Diffusion function.
        if flag,
            cW = 1./(1 + (nablaW/kappa).^2);
            cE = 1./(1 + (nablaE/kappa).^2);
        else
            cW = exp(-(nablaW/kappa).^2);
            cE = exp(-(nablaE/kappa).^2); 
        end;

        % Discrete PDE solution.
        diff_sig = diff_sig + ...
                   delta_t*(...
                   (1/(dx^2))*cW.*nablaW + (1/(dx^2))*cE.*nablaE);
        
        % Iteration warning.
%         fprintf('\rIteration %d\n',t);
end