%%%%% Douglas Rachford in 2 ways
% DR splitting solves problems like
%     min f(x) + g(x)
%     x

A = rand( 7, 3 ) * 100;
x = rand( 3, 1 ) * 100;
b = A * x;
lambda = 15;
f = @(x) lambda * norm( x, 1 );
g = @(x) 0.5 * norm( A*x - b, 2 )^2;
proxf = @(x,t) softThresh( x, lambda * t );
proxg = @(y,t) proxL2Sq( y, 1, b, A );  % TODO:  should there be a t input here?
t = 1d-3;  % admm step size

%%% method 1
% the standard DR iteration
% xk = proxf( zk, gamma )
% yk = proxg( 2xk - zk, gamma )
% zk+1 = zk + 2alpha(yk - xk)
zk = zeros( size( x ) );

gamma = 1;
alpha = 0.5;
for i = 1:15
  xk = proxf(zk, gamma);
  yk = proxg(2*xk - zk, gamma);
  zkp1 = zk + 2*alpha*(yk - xk);
  zk = zkp1;
  fprintf('iter %d f(x) = %f g(x) = %f \n', i, f(zk), g(zk));
end

norm(x - zk) / norm(x)
f(zk) + g(zk)
f(x) + g(x)

%%% method 2
% writing DR as a fixed point
% zk+1 = zk + alpha( R_gamma,g [ R_gamma,f zk] - zk )
% where R_gamma,g (x)  = 2proxg(x, gamma) - x

zk = zeros( size( x ) );
gamma = 1;
alpha = 0.5;

Rg = @(x, gamma) 2 * proxg(x, gamma) - x;
Rf = @(x, gamma) 2 * proxf(x, gamma) - x;

for i = 1:15
  zkp1 = zk + alpha* (Rg( Rf( zk, gamma ), gamma) - zk);
  zk = zkp1;
  fprintf('iter %d f(x) = %f g(x) = %f \n', i, f(zk), g(zk));
end
norm(x - zk) / norm(x)
f(zk) + g(zk)
f(x) + g(x)

%%% linesearch
% use DR as fixed point and use linesearch from paper

gamma = 1/3;
eps = 0.003;
alpha_change = 1/1.4;
S = @(x) Rg( Rf( x, gamma ), gamma );

xk = zeros( size( x ) );
alpha_bar = alpha;

for i = 1:15
  rk = S(xk) - xk;
  xk_bar = xk + alpha_bar*rk;
  rk_bar = S(xk_bar) - xk_bar;

  alpha_k = 50;
  subiter = 0;
  while true
    subiter = subiter + 1;
    fprintf('subiter %d\n', subiter)
    xkp1 = xk + alpha_k*rk;
    rkp1 = S(xkp1) - xkp1;
    if norm(rkp1) < (1-eps)*norm(rk_bar)
      xk = xkp1;
      break
    end
    alpha_k = alpha_k*alpha_change;
    if alpha_k < alpha_bar
      alpha_k = alpha_bar;
      xkp1 = xk + alpha_k*rk;
      xk = xkp1;
      break
    end
  end
end

norm(x - xk) / norm(x)
f(xk) + g(xk)
f(x) + g(x)