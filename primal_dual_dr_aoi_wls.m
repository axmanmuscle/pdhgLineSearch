function [xStar, iters, alphas, objVals] = primal_dual_dr_aoi_wls(x0,proxf,proxgconj, f, g, maxIter)

%%% parameters
alpha_bar = 0.5; % alpha_bar
gamma = 3; % gamma for prox operators
eps = 0.03; % eps for (1 - eps) || rbar_k || in linesearch
tol = 1e-7; % tolerance for exit criterion
alpha_change = 1/1.4; % factor for change in alpha during linesearch

  function out = Rf( phi )
    out = 2 * proxf( phi, gamma ) - phi;
  end

  function out = Rg2( phi )
    out = 2 * proxgconj( phi, 1/gamma ) - phi;
  end

  function out = S ( in )
    out = -gamma * Rg2( (1/gamma) * Rf(in) );
  end

%Rf = @(phi) 2*proxf(phi, gamma) - phi;
% Rgconj = @(phi) 2*proxgconj(phi, gamma) - phi;
%Rg2 = @(phi) 2*proxgconj(phi, 1/gamma) - phi;

% S = @(phi) Rgconj(Rf(phi));
%S = @(phi) -gamma * Rg2( (1/gamma) * Rf(phi) );

iters = zeros([maxIter size(x0)]);
alphas = zeros(maxIter, 1);
objVals = zeros([maxIter, 1]);

xk = x0;

for i = 1:maxIter
  rk = S(xk) - xk;
  xk_bar = xk + alpha_bar * rk;
  rk_bar = S(xk_bar) - xk_bar;

  alpha_k = 50;
  subiter = 0;
  while true
    subiter = subiter + 1;
%     fprintf('subiter %d\n', subiter);
    xkp1 = xk + alpha_k * rk;
    rkp1 = S(xkp1) - xkp1;
    if norm(rkp1) < (1-eps) * norm(rk_bar)
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

  iters(i, :) = xk;

  objVals(i) = f(xk) + g(xk);

  %if mod(i, 3) == 1
    fprintf('iter: %d   objVal: %d   res: %d  alpha: %d\n', i, objVals(i), norm(rkp1(:)), alpha_k );
  %end

end
xk = proxf(xk, gamma);
xStar = xk;
end