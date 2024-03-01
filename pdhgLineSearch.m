%%%%%%%%%%%
% making a line search for our PDHG method for compressed sensing
% reconstruction with homodyne detection. we'll use the paper to go from
% PDHG -> D-R, then use the averaged operator iteration form of D-R to code
% up our line search

A = randn(10, 30);
x = rand(30, 1);

b = A*x;

norm(A\b - x);

%%%%%%%%%%
% let's use the following problem
%           min || x ||_1 + || Ax - b ||_2^2
%            x
% we'll solve this with PDHG and then with DR?