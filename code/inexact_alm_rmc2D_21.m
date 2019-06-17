function [X_hat, E_hat, O_hat, iter] = inexact_alm_rmc2D_21(D, Sigma_bar, lambda, tol, maxIter)
% Jun 2019
% This matlabcode implements the robust matrix completion with column-sparse
% gross corruption.
% 
% [X_hat, E_hat, O_hat, iter] = inexact_alm_rmc2D_21(D, Sigma_bar, lambda, tol, maxIter)
% returns the estimated complete low rank matrix X_hat, and fiber-sparse
% corruption matrix E_hat. -1 * O_hat is estimation for what should have been
% in D for missing entries. iter is the number of total iterations.
%
% D - observation data in matrix format. The outlier fiber is assumed to be
% arranged along the columns. 
%
% Sigma_bar - index for unobserved entries. It is a matrix the same size as
% D, with value one corresponding to missing entries in D, and zero
% otherwise.
%
% lambda - weight on fiber-sparse corruption matrix in the cost function.
%        - Larger lambda results in sparser corruption matrix.
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
%
% maxIter - maximum number of iterations
%         - DEFAULT 1000, if omitted or -1.
% 
% solves the optimization problem:
% min_{X, E} |X_(i)|_* + lamnda * |E|_{2,1}
% s.t        D = X + E + O;
%            O .* ~Sigma_bar = 0.
% 
% Yue Hu, Jun 2019. Questions? yue.hu@vadnerbilt.edu;
% Daniel B. Work (dan.work@vanderbilt.edu)
%
%
% This code is inspired by the work of Zhouchen Lin, Risheng Liu, and 
% Zhixun Su, Linearized Alternating Direction Method with Adaptive Penalty 
% for Low Rank Representation, NIPS 2011.
%
% We are grateful for the permission of Professor John Wright to modify,
% extend and share their code of Augmented Lagrange Multiplier (ALM) Method
% for robust matrix PCA, for the purpose of academic research. Their
% original code, copyrighted by Perception and Decision Lab, University of
% Illinois at Urbana-Champaign, and Microsoft Research Asia, Beijing, can
% be found at
% https://people.eecs.berkeley.edu/~yima/matrix-rank/sample_code.html.
%


addpath PROPACK;

[m n] = size(D);

if nargin < 2
    lambda = 1 / sqrt(m);
end

if nargin < 3
    tol = 1e-7;
elseif tol == -1
    tol = 1e-7;
end

if nargin < 4
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

% initialize
Y = D;
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

% Sigma = D~= 0;
% Sigma_bar = D==0;  % projection matrix for non-observed entries
X_hat = zeros( m, n);
E_hat = zeros( m, n);
O_hat = zeros( m, n);

mu = 1.25/norm_two; % this one can be tuned
mu_bar = mu * 1e7;
rho = 1.5;          % this one can be tuned
d_norm = norm(D, 'fro');

iter = 0;
total_svd = 0;
converged = false;
stopCriterion = 1;
sv = 10;

while ~converged       
    iter = iter + 1;
    
    % update E
    temp_T = D - X_hat + (1/mu)*Y - O_hat;
    for j = 1:size(temp_T,2)
        E_hat(:,j) = temp_T(:,j) * max(0,1-lambda/norm(temp_T(:,j)));
    end    
%     E_hat = E_hat .* Sigma
%     col_X = ~any(E_hat);  %find all nonzero coloumns of E
%     Col_X = repmat(col_X,m,1); 

    % update X
    if choosvd(n, sv) == 1
        [U S V] = lansvd(D - E_hat + (1/mu)*Y - O_hat, sv, 'L');
    else
        [U S V] = svd(D - E_hat + (1/mu)*Y - O_hat, 'econ');
    end
    diagS = diag(S);
    svp = length(find(diagS > 1/mu));
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end
    
    X_hat = U(:, 1:svp) * diag(diagS(1:svp) - 1/mu) * V(:, 1:svp)';    
%     X_hat = X_hat .* Col_X;
    
    total_svd = total_svd + 1;
    
        
    % update O
    O_temp = D - X_hat - E_hat + (1/mu)*Y;
    O_hat = Sigma_bar .* O_temp;
    
    
    % update Y
    Z = D - X_hat - E_hat - O_hat;
    
    Y = Y + mu*Z;
    mu = min(mu*rho, mu_bar);
        
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / d_norm;
    if stopCriterion < tol
        converged = true;
    end    
    
    if mod( total_svd, 10) == 0
%         disp(['#svd ' num2str(total_svd) ' r(A) ' num2str(rank(X_hat))...
%             ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
%             ' stopCriterion ' num2str(stopCriterion)]);
    end    
    
    if ~converged && iter >= maxIter
%         disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end

% Treats the entire coloumn as outlier once it is corrupted. 
% i.e. set the coloumn of low rank matrix to zero if the corresponding
% Sparse matrix have value.
%
% E_hat = E_hat .* Sigma
col_X = ~any(E_hat);  %find all nonzero coloumns of E
Col_X = repmat(col_X,m,1); 

X_hat = X_hat .* Col_X;
E_hat = D .* (~Col_X);


end
