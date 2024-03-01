
function test_pdhg_aoi()
  
%   if numel( A ) == 0
%     applyA = @(x) x;
%     applyAT = @(x) x;
%   elseif isnumeric( A )
%     applyA = @(x) A * x;
%     applyAT = @(y) A' * y;
%   else
%     applyA = @(x) A( x, 'notransp' );
%     applyAT = @(x) A( x, 'transp' );
%   end

  alpha0 = 0.95;
  gamma = 2;

  rng(20240124)

  %%% data
  A = randn(20, 20);
  x = randn(20, 1);
  b = A*x;
  
  %%% functions
  f = @(x) norm(A*x - b, 2)^2;
  % proxf = @(y,t) proxL2Sq( y, 1, b, A );
  applyA = @(x) A * x;
  applyAT = @(y) A' * y;
  
  P = inv(A'*A + eye(size(A))); % compute inverse once and use thousands of times
  
  proxf = @(y, t) P * (y - b);
  proxgConj = @(x, t) x - t*proxg(x/t, 1/t);
  
  Rg = @(x, gamma) 2 * proxg(x, gamma) - x;
  Rf = @(x, gamma) 2 * proxf(x, gamma) - x;

  alphak = alpha0;

  %%% Averaged operator iteration is of the form
  % z_k+1 = (1 - alpha) z_k + alpha S z_k

  theta = 2;
  gamma = .25;
  alpha = 0.05;
 
  n = size(x, 1);
  x0 = zeros(n,1);

  phik = x0;
  pk = zeros(size(x));
  zk = 0;
  optIter = 1;
  N = 4000;
  while optIter < N
    optIter = optIter + 1;

    ak = phik(1:n);

    akp1 = proxf(ak - pk, gamma);

    Sphik = (theta / gamma)*applyA(2*akp1 - ak) + zk;
    Sphik = Sphik - 2*proxgConj(Sphik, theta/gamma);
    Sphik = Sphik * alpha;

    zkp1 = (1 - alpha)*zk + (theta / gamma)*applyA(akp1); %???
    zkp1 = zkp1 - (1 - alpha)*(theta/gamma)*applyA(ak) - alpha*Sphik;

    phikp1 = (1 - alpha)*phik + alpha*Sphik;

    pkp1 = phikp1(1:n) - akp1;


    phik = phikp1;
    zk = zkp1;
    pk = pkp1;
    
  end

  if nargout > 1, objValues = objValues( 1 : optIter ); end
  if nargout > 2, relDiffs = relDiffs( 1 : optIter ); end

  xStar = x;
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
