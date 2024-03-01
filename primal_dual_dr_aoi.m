function [xStar, iters] = primal_dual_dr_aoi(x0,proxf,proxgconj)

Rf = @(phi) 2*proxf(phi, 1) - phi;
Rg = @(phi) phi - 2*proxgconj(phi, 1);
Rgconj = @(phi) 2*proxgconj(phi, 1) - phi;

S = @(phi) -1*Rgconj(Rf(phi));
Sstar = @(phi) Rgconj(Rf(phi));

maxIter = 100;
iters = zeros([maxIter size(x0)]);

xk = x0;
alpha = 0.5;
for i = 1:maxIter
%   xkp1 = (1-alpha)*xk + alpha*S(xk);
%   xkp1 = 0.5*(S(xk) + xk);
  xkp1 = 0.5*xk - 0.5*Sstar(xk);
  xk = proxf(xkp1);

  iters(i, :) = proxf(xk, 1);

end

xStar = xk;
end