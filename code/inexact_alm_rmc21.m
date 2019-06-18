function [X_hat, E_hat, O_hat,iter] = inexact_alm_rmc21(D,Sigma_bar, lambda, tol, maxIter)
% Jun 2019
% This matlabcode implements the robust tensor completion with fiber-sparse
% gross corruption.
% 
% [X_hat, E_hat, O_hat,iter] = inexact_alm_rmc21(D,Sigma_bar, lambda, tol, maxIter)
% returns the estimated complete low rank tensor X_hat, and fiber-sparse
% corruption tensor E_hat. -1 * O_hat is estimation for what should have been
% in D for missing entries. iter is the number of total iterations.
%
% D - observation data in tensor format. The outlier fiber is assumed to be
% arranged along the first mode. 
%
% Sigma_bar - index for unobserved entries. It is a tensor the same size as
% D, with value one corresponding to missing entries in D, and zero
% otherwise.
%
% lambda - weight on fiber-sparse corruption tensor in the cost function.
%        - Larger lambda results in sparser corruption tensor.
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
%
% maxIter - maximum number of iterations
%         - DEFAULT 1000, if omitted or -1.
% 
% solves the optimization problem:
% min_{X, E} sum_i(|X_(i)|_*) + lamnda * |E_(1)|_{2,1}
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
addpath PROPACK/tensor_toolbox-master ;

if nargin < 4
    tol = 1e-7;
elseif tol == -1
    tol = 1e-7;
end

if nargin < 5
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

% initialize
Y0 = D;
Y_mat = tenmat(Y0,1);
norm_two = lansvd(Y_mat.data, 1, 'L');
norm_inf = norm( Y_mat.data(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y0 = Y0 / dual_norm;

D_mode = ndims(D); 
Y = cell(D_mode,1);
for i = 1:D_mode
    Y{i} = Y0;
end


X_hat = tenzeros( size(D));
E_hat = tenzeros( size(D));
O_hat = tenzeros( size(D));
mu = 1/norm_two; % this one can be tuned
mu_bar = mu * 1e7;
rho = 1.5;          % this one can be tuned
d_norm = norm(D);

Sigma_bar = tensor(Sigma_bar);

iter = 0;
total_svd = 0;
converged = false;
stopCriterion = 1;


D_size = size(D);
sv0 = min(D_size);
sv = cell(D_mode,1);   %number of singular values desired 
for i = 1:D_mode
    sv{i} = sv0;
end

X_all = cell(D_mode,1);


while ~converged       
    iter = iter + 1;
%     fprintf('mu is %.2f\n',mu)
    % update E
    Y_mean = tenfun(@sum,Y{:})./D_mode;
    temp_T = D - X_hat + (1/mu)*Y_mean - O_hat;
    temp_Tm = tenmat(temp_T,1);
    for j = 1:size(temp_Tm,2)
        temp_Tm(:,j) = temp_Tm(:,j) * max(0,1-lambda/(D_mode*mu*norm(temp_Tm(:,j))));
%         E_m(:,j) = temp_Tm(:,j) * max(0,1-lambda/norm(temp_Tm(:,j)));
    end 
    E_hat = tensor(temp_Tm);
%     E_hat = E_hat .* ~Sigma_bar;
    
    % update X(i)
    for i = 1:D_mode
        temp = D - E_hat + (1/mu)*Y{i} - O_hat;
        temp_mat = tenmat(temp,i);
        n_ = size(temp_mat);
        n = min(n_);
        sv{i} = min(sv{i},n);
        
        if choosvd(n, sv{i}) == 1            
            [U, S, V] = lansvd(temp_mat.data, sv{i}, 'L');
        else
            [U, S, V] = svd(temp_mat.data, 'econ');
        end      

        diagS = diag(S);
        svp = length(find(diagS > 1/mu));
        if svp < sv{i}
            sv{i} = min(svp + 1, n);
        else
            sv{i} = min(svp + round(0.05*n), n);
        end

        temp_mat(:,:) = U(:, 1:svp) * diag(diagS(1:svp) - 1/mu) * V(:, 1:svp)'; 
    
        X_all{i} = tensor(temp_mat);
    end
    
    X_hat = tenfun(@sum,X_all{:})./D_mode;  
        
    
    total_svd = total_svd + 1;
    
    % update O
    O_temp = D - X_hat - E_hat + (1/mu)*Y_mean;
    O_hat = Sigma_bar .* O_temp;
   
    % update Y
    for i = 1:D_mode
        Z = D - X_all{i} - E_hat - O_hat;
        Y{i} = Y{i} + mu*Z;
    end
    
    mu = min(mu*rho, mu_bar);
        
    %% stop Criterion    
    stopCriterion = norm(Z) / d_norm;
    if stopCriterion < tol
        converged = true;
    end    
   
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end

% Treats the entire coloumn as outlier once it is corrupted. 
% i.e. set the coloumn of low rank matrix to zero if the corresponding
% Sparse matrix have value.

E_m = tenmat(E_hat,1);
col_X = any(double(E_m) ~= 0);  %find index all nonzero coloumns of E

X_m = tenmat(X_hat,1);
E_m = tenmat(D,1);
X_m(:,col_X) = 0;
E_m(:,~col_X) = 0;

X_hat = tensor(X_m);
E_hat = tensor(E_m);





end
