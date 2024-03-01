function [xStar, iters] = primal_dual_dr_3(a0,proxf,proxgconj, A, B, theta)

if isnumeric(A)
  applyA = @(x) A*x;
  applyAt = @(x) A'*x;
  applyB = @(x) B*x;
  applyBt = @(x) B'*x;
else
  applyA = @(x) A(x, 'notransp');
  applyAt = @(x) A(x, 'transp');
  applyB = @(x) B(x, 'notransp');
  applyBt = @(x) B(x, 'transp');
end

maxIter = 100;
iters = zeros([maxIter size(a0)]);

n = size(a0, 1);
ak = a0;
ck = zeros(size(ak));
dk = zeros(size(ak));

for i = 1:maxIter
  akp1 = proxf(ck, 1);
  ckp1 = akp1 - applyAt(proxgconj(theta*applyA(2*akp1 - ck), 1));

  tmp1 = proxgconj(theta*applyA(2*akp1 - ck) - applyB(dk), 1);
  tmp2 = applyAt(tmp1);
  tmp3 = applyBt(tmp1);

  update_actual = [akp1;dk] - [tmp2; tmp3];
  ckp1_actual = update_actual(1:n);
  dkp1 = update_actual(n+1:end);

%   if norm(ckp1_actual - ckp1) > 1e-4
%     error('Problem with your assumptions in the PD-DR-3 method')
%   end

  ak = akp1;
  ck = ckp1;
  dk = dkp1;

  iters(i, :) = ak;

end

xStar = ak;
end