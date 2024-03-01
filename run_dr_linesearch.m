%%% linesearch test
% use DR as fixed point and use linesearch from paper
function [] = run_dr_linesearch()

rng(20240124)
%%% data
A = randn(20, 20);
x = randn(20, 1);
b = A*x;

%%% functions
f = @(x) norm(A*x - b, 2)^2;
% proxf = @(y,t) proxL2Sq( y, 1, b, A );

P = inv(A'*A + eye(size(A))); % compute inverse once and use thousands of times

proxf = @(y, t) P * (y - b);

Rg = @(x, gamma) 2 * proxg(x, gamma) - x;
Rf = @(x, gamma) 2 * proxf(x, gamma) - x;

%%% parameters
alpha = 0.5;
gamma = 3;
eps = 0.03;
tol = 1e-3;
alpha_change = 1/1.4;
max_iter = 5000;
S = @(x) Rg( Rf( x, gamma ), gamma );
iters = zeros(max_iter, 1);
alphas = zeros(max_iter, 1);

xk = zeros( size( x ) );
alpha_bar = alpha;

for i = 1:max_iter
%   if mod(i, 10) == 0
%     fprintf('iter %d\n', i)
%   end
  rk = S(xk) - xk;
  xk_bar = xk + alpha_bar*rk;
  rk_bar = S(xk_bar) - xk_bar;

  iters(i) = norm(rk);
  if norm(rk) < tol
    disp('break')
    break
  end

  alpha_k = 50;
  subiter = 0;
  while true
    subiter = subiter + 1;
    %fprintf('subiter %d\n', subiter);
    xkp1 = xk + alpha_k*rk;
    rkp1 = S(xkp1) - xkp1;
    if norm(rkp1) < (1-eps)*norm(rk_bar)
      xk = xkp1;
      alphas(i) = alpha_k;
      break
    end
    alpha_k = alpha_k*alpha_change;
    if alpha_k < alpha_bar
      alpha_k = alpha_bar;
      xkp1 = xk + alpha_k*rk;
      xk = xkp1;
      alphas(i) = alpha_k;
      break
    end
  end
end

[xstar, iters1, alphas1] = dr_wLS(zeros(size(x)), proxf, @proxg, 5000);

norm(xstar - xk)

norm(iters - iters1, 'fro');

norm(S(xk) - xk)
semilogy(iters)
norm(x - xk) / norm(x)
f(xk) + g(xk)
f(x) + g(x)
end

function gx = g(x)
  if sum(x < 0) > 0
    gx = Inf;
  else
    gx = 0;
  end
end

function out = proxg(x, t)
  out = x;
  out(x < 0) = 0;
end