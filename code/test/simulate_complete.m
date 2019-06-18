function [D, L ,S] = simulate_complete(D1,D2,D3,C1,C2,C3,ratio_s )
% Jun 2019
% [D, L ,S] = simulate_complete(D1,D2,D3,C1,C2,C3,ratio_s ) creates tensor
% with deminsion(D1,D2,D3), Tucker rank (C1,C2,C3) that is fiber-wise corrupted.
% D is the observation tensor, L is the ground truth low rank tebsor, and S
% is the fiber-sparse tensor.
%
% ratio_s - corruption ratio, default 0.05.
% D - The simulated obsevsation data
% L - The low rank matrix
% S - The Sparse matrix
% 
% Yue Hu, Jun 2019. Questions? yue.hu@vadnerbilt.edu;

addpath PROPACK;
addpath PROPACK/tensor_toolbox-master
% default: 5% corruoted column
if nargin < 7
    ratio_s = 0.05;
end

% LOW RANK TENSOR
C = tensor(randn(C1,C2,C3));
A1 = RandOrthMat(D1,C1);
A2 = RandOrthMat(D2,C2);
A3 = RandOrthMat(D3,C3);

L = ttm(C,{A1,A2,A3});

% % check the actuall rank
% T = hosvd(L,1e-6);
% coresize = size(T.core);   % Check size of core
% fprintf('the tucker rank for low rank tensor L is(%d,%d,%d)\n',coresize(1),coresize(2),coresize(3))


% sparse tensor, column sparse
S_ = rand(size(L))*1;  %? grossly corrsupted
Col_ = rand(1,D2*D3);
Col_S = Col_<= ratio_s; % ratio_s corrupted column
Col_S = repmat(Col_S,D1,1);
Col_S = tensor(Col_S,size(S_));
S = S_ .* Col_S;

% Corresponding column for L is 0
All = tenones(size(L));
Col_L = All - Col_S;
L = L .* Col_L;

% Observation matrix
D = S + L;

end