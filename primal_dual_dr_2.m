function [xStar, iters] = primal_dual_dr_2(x0,proxf,proxgconj)

maxIter = 100;
iters = zeros([maxIter size(x0)]);

xk = x0;
yk = zeros(size(xk));
for i = 1:maxIter
  xkp1 = proxf(yk, 1);
  ykp1 = xkp1 - proxgconj(2*xkp1 - yk, 1); 

  xk = xkp1;
  yk = ykp1;

  iters(i, :) = xk;

end

xStar = xk;
end