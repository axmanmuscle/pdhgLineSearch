function run_test_prob_dr()
b = [5;5;5];
eps = 0.3;
n = size(b, 1);

%% problem is 
% min f(x) + g(Ax)
%A = [2 0 1;1 2 3; 8 4 6];
A = randi(15, [9 3]);
theta = 1/norm(A)^2 - 1e-6;
% theta = 2;
tau = 0.1;
Bt = chol((1/theta)*eye(9) - A*A');
B = Bt';

sigma = theta/tau;

g = @(x) norm(x, 1);
f = @(x) 0;

gtilde = @(x) norm(A*x(1:n) + B*x(n+1:end), 1);
ftilde = @(x) full_f(x(1:n), b, eps);
pf = @(x, t) [proxF(x(1:n), b, eps); zeros([9 1])];
pgstar = @(x, t) proxConjL1(x, 1, t);
pgtilde = @(x, t) x - tau * [A'; B'] * pgstar(sigma*(A*x(1:n) + B*x(n+1:end)), sigma);
pgtildestar = @(x, t) x - pgtilde(x, t);

pg3 = @(x, t) tau * [A'; B'] * pgstar(sigma*(A*x(1:n) + B*x(n+1:end)), sigma);

% x0 = [-20;-10; -30; 0; 0; 0];
x0 = [-20;-10;-30;zeros([9 1])];

% xStar = douglasRachford(x0, pf, pgtilde, 1);
[xStar2, driters] = primal_dual_dr(x0, pf, pgtildestar);
[xStar3, dr2iters] = primal_dual_dr_2(x0, pf, pgtildestar);
[xStar4, aoiiters] = primal_dual_dr_aoi(x0, pf, pgtildestar);
[xStar6, aoiiters2, alphas] = primal_dual_dr_aoi_wls(x0, pf, pg3, f, g, 300);

pfnew = @(x, t) proxF(x(1:n), b, eps);
[xStarpdhg, objvals, reldiffs] = pdhg(x0(1:3), pfnew, pgstar, tau, 'A', A, ...
  'f', f, 'g', g, 'normA', norm(A), 'verbose', true, 'tol', 1e-10);
[xStar5, dr3iters] = primal_dual_dr_3(x0(1:n), pfnew, pgstar, A, B, theta);


xStar

end

function out = full_f(x, b, eps)
  if norm(x - b, 2) <= eps
    out = 0;
  else
    out = Inf;
  end
end

function out = proxF(x, b, eps)
  if norm(x - b, 2) < eps
    out = x;
  else
    xperp = (x - b)/norm(x-b,2) * eps;
    out = b + xperp;
  end
end