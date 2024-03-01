function [xStar, iters] = primal_dual_dr(x0,proxf,proxgconj)

maxIter = 100;

xk = x0;
zk = zeros(size(xk));

iters = zeros([maxIter size(xk)]);
for i = 1:maxIter
  xkp1 = proxf(xk - zk, 1);
  zkp1 = proxgconj(zk + 2*xkp1 - xk, 1); 

  xk = xkp1;
  zk = zkp1;

  iters(i, :) = xk;

end

xStar = xk;
end