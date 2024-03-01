%%% adaptive DR

function run_adaptive_DR()

%%% data
A = randn(50, 50);
x = randn(50, 1);
b = A*x;

%%% functions
f = @(x) norm(A*x - b, 2)^2;
proxf = @(y,t) proxL2Sq( y, 1, b, A );

% P = inv(A'*A + eye(size(A))); % compute inverse once and use thousands of times
% proxf = @(y, t) P * (y - b);

x0 = randn(size(x));
t0 = 1;

xn = x0;
tn = t0;

max_iters = 300;
iters = zeros(max_iters, 1);
for i = 1:max_iters
  un = proxg(xn, tn);
  kn = norm(un) / norm(un - xn);

  tn = kn*tn;

  vn = proxf((1+kn)*un - kn*xn, tn);

  xn = vn + kn*(xn - un);

  iters(i) = f(proxg(xn, tn));
end

figure; plot(iters);

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