%%% PDHG
% minimize 0.5*|| Ax - b ||^2

function run_pdhg()

%%% data
rng(20231130)
A = randn(500, 500);
x = randn(500, 1);
b = A*x;

L = norm(A);
tau = 1/L;
sigma = 1/L;
theta = 1;
gamma = 1/L;

%%% functions
f = @(x) 0.5*norm(A*x - b, 2)^2;
proxf = @(y,t) proxL2Sq( y, 1, b, A );

% P = inv(A'*A + eye(size(A))); % compute inverse once and use thousands of times
% proxf = @(y, t) P * (y - b);

x0 = zeros(size(x));
p0 = zeros(size(x));

xbar = x0;

xn = x0;
pn = p0;

max_iters = 1000;
iters = zeros(max_iters, 1);

for i = 1:max_iters
  pn = (pn + sigma*(A*xbar - b)) / (1 + sigma);

  x_old = xn;

  xn = xn - tau*A'*pn;

%   %%% update params
  theta = 1/sqrt(1+2*gamma*tau);
  tau = theta*tau;
  sigma = sigma / theta;
  
  %%% update dual variable
  xbar = xn + theta*(xn - x_old);

  iters(i) = f(xbar);
  fprintf('iter %d error %f\n', i, f(xbar))
end

figure; semilogy(iters);

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