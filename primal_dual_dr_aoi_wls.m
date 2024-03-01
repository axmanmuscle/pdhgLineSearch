function [xStar, iters, alphas] = primal_dual_dr_aoi_wls(x0,proxf,proxgconj)

%%% parameters
alpha_bar = 0.5; % alpha_bar
gamma = 0.1; % gamma for prox operators
eps = 0.03; % eps for (1 - eps) || rbar_k || in linesearch
tol = 1e-7; % tolerance for exit criterion
alpha_change = 1/1.4; % factor for change in alpha during linesearch

Rf = @(phi) 2*proxf(phi, gamma) - phi;
Rgconj = @(phi) 2*proxgconj(phi, gamma) - phi;
Rg2 = @(phi) phi - 2*proxgconj(phi, gamma);

% S = @(phi) Rgconj(Rf(phi));
S = @(phi) Rg2(Rf(phi));

maxIter = 100;
iters = zeros([maxIter size(x0)]);
alphas = zeros(maxIter, 1);

xk = x0;

for i = 1:maxIter
  if mod(i, 10) == 1
    fprintf('iter %d\n', i)
  end
  rk = S(xk) - xk;
  xk_bar = xk + alpha_bar*rk;
  rk_bar = S(xk_bar) - xk_bar;

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

  xk = proxf(xk, gamma);

  iters(i, :) = xk;

end

xStar = xk;
end