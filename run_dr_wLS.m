%%% linesearch test
% use DR as fixed point and use linesearch from paper
function [] = run_dr_wLS()

%%% data
A = randn(1000, 1000);
x = randn(1000, 1);
b = A*x;

%%% solve min || Ax - b ||_2^2 + lambda || x ||_1
lambda = 5;

%%% functions
f = @(x) norm(A*x - b, 2)^2;
g = @(x) lambda*norm(x, 1);
% proxf = @(y,t) proxL2Sq( y, 1, b, A );
proxg = @(x, t) softThresh(x, lambda*t);

P = inv(A'*A + eye(size(A))); % compute inverse once and use thousands of times
proxf = @(y, t) P * (y - b);

[xstar, iters1, alphas1] = dr_wLS(zeros(size(x)), proxf, proxg, 4000);

% disp(xstar);
f(xstar) + g(xstar)

figure; semilogy(iters1);
figure; plot(alphas1);


end

% function gx = g(x, lambda)
%   if sum(x < 0) > 0
%     gx = Inf;
%   else
%     gx = 0;
%   end
% end
% 
% function out = proxg(x, t)
%   out = x;
%   out(x < 0) = 0;
% end