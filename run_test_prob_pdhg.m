function run_test_prob_pdhg()
b = [5;5];
eps = 0.3;

%% problem is 
% min f(x) + g(Ax)
A = [2 0; 0 3];

g = @(x) norm(x, 1);
f = @(x) full_f(x, b, eps);
pf = @(x, t) proxF(x, b, eps);
pgstar = @(x, t) proxConjL1(x, 1, t);

x0 = [7;3];
tau = 0.1;
[xStar,objValues,relDiffs] = pdhg(x0, pf, pgstar, tau, 'A', A, 'normA', norm(A), 'f', f, 'g', g);

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