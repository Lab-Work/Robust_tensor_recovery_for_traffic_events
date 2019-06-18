function [D, L ,S] = simulate_4D(D1,D2,D3,D4,C1,C2,C3,C4,ratio_s )
% [D, L ,S] = simulate3(D1,D2,D3,D4,C1,C2,C3,C4,ratio_s) creates tensor
% with 4D deminsion(D1,D2,D3,D4), Tucker rand(C1,C2,C3,C4) 
% that is fiber-wise corrupted for ratio_s of the fibers.
% D is the observation tensor, L is the ground truth low rank tebsor, and S
% is the fiber-sparse tensor.
%
% ratio_s - corruption ratio, default 0.05.
% D - The simulated obsevsation data
% L - The low rank matrix
% S - The Sparse matrix
% 

addpath PROPACK;
addpath PROPACK/tensor_toolbox-master
% default: 5% corruoted column
if nargin < 9
    ratio_s = 0.05;
end

% LOW RANK TENSOR
C = tensor(randn(C1,C2,C3,C4));
A1 = RandOrthMat(D1,C1);
A2 = RandOrthMat(D2,C2);
A3 = RandOrthMat(D3,C3);
A4 = RandOrthMat(D4,C4);

L = ttm(C,{A1,A2,A3,A4});

% sparse tensor, column sparse
S_ = rand(size(L))*1;  %? grossly corrsupted
Col_ = rand(1,D2*D3*D4);
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