
function [xStar,objValues,relDiffs] = pdhg_aoi( phi, proxf, proxgConj, tau, varargin )
  % [xStar,objValues,relDiffs] = pdhg( x, proxf, proxgConj, tau [, ...
  %   'A', A, 'f', f, 'g', g, 'N', N, 'normA', normA, 'sigma', sigma, ...
  %    'lambda', lambda, 'printEvery', printEvery, 'theta', theta, ...
  %    'tol', tol, 'verbose', verbose, 'z', z ] )
  %
  % Implements the Primal Dual Hybrid Gradient (Chambolle-Pock) method that
  % solves problems of the form:  minimize f( x ) + g( A x ), in the form
  % of an averaged operator iteration ustilizing Boyd's line search.
  %
  % Inputs:
  % x - initial guess
  % proxf - a function handle for the proximal operator of
  % proxgConj - a function handle for the proximal operator of the conjugate of g
  %
  % Optional Inputs:
  % A - if A is not provided, it is assumed to be the identity
  % f - to determine the objective values, f must be provided
  % g - to determine the objective values, g must be provided
  % N - the number of iterations that ADMM will perform (default is 100)
  % normA - the matrix induced 2-norm of A.  Could be determined with norm or, for
  %   large sparse matries, estimated with normest or powerIteration.
  % lambda - relaxation parameter
  % theta - acceleration parameter
  % verbose - true or false
  % z - initial value of dual variable
  %
  % Outputs:
  % xStar - the optimal point
  %
  % Optional Outputs:
  % objValues - a 1D array containing the objective value of each iteration
  % relDiffs - a 1D array containing the relative difference of each iteration
  %
  % Written by Nicholas Dwork - Copyright 2019
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.

  if nargin < 1
    disp( 'Usage;  ' );
    disp( '  [xStar,objValues] = pdhg( x, proxf, proxgConj, tau [, ... ');
    disp( '  ''A'', A, ''f'', f, ''g'', g, ''N'', N, ''normA'', normA, ''sigma'', sigma, ...' );
    disp( '  ''lambda'', lambda, ''printEvery'', printEvery, ''theta'', theta, ...' );
    disp( '  ''tol'', tol, ''verbose'', verbose, ''z'', z ] )' );
    if nargout > 0, xStar = []; end
    if nargout > 1, objValues = []; end
    return;
  end

  defaultTol = 1d-4;

  p = inputParser;
  p.addParameter( 'A', [] );
  p.addParameter( 'f', [] );
  p.addParameter( 'g', [] );
  p.addParameter( 'N', 100, @ispositive );
  p.addParameter( 'normA', [], @(x) x>0 );
  p.addParameter( 'sigma', [], @ispositive );
  p.addParameter( 'theta', 1, @(x) x >= 0 && x <= 1 );
  p.addParameter( 'tol', defaultTol, @(x) numel(x) == 0 || ispositive(x) );
  p.addParameter( 'lambda', 1, @(x) x >= 0 && x <= 2 );
  p.addParameter( 'verbose', false, @(x) islogical(x) || x==1 || x==0 );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'z', [], @isnumeric );
  p.parse( varargin{:} );
  A = p.Results.A;
  f = p.Results.f;
  g = p.Results.g;
  N = p.Results.N;
  normA = p.Results.normA;
  sigma = p.Results.sigma;
  theta = p.Results.theta;
  tol = p.Results.tol;
  lambda = p.Results.lambda;
  printEvery = p.Results.printEvery;
  verbose = p.Results.verbose;
  z = p.Results.z;
  
  if numel( A ) == 0
    applyA = @(x) x;
    applyAT = @(x) x;
  elseif isnumeric( A )
    applyA = @(x) A * x;
    applyAT = @(y) A' * y;
  else
    applyA = @(x) A( x, 'notransp' );
    applyAT = @(x) A( x, 'transp' );
  end

  if nargout > 1,  objValues = zeros( N, 1 ); end
  if nargout > 2,  relDiffs = zeros( N, 1 ); end

  if numel( sigma ) == 0  && numel( A ) > 0
    if numel( normA ) == 0
      error( 'If an A is supplied, you must supply sigma or normA' );
    end
    if normA == 0, return; end
    sigma = ( 0.99 / normA^2 ) / tau;
  end

  calculateRelDiffs = false;
  if nargout > 2  ||  ( numel( tol ) > 0  &&  tol < Inf )
    calculateRelDiffs = true;
  end

  maxAbsX = max( abs( x(:) ) );
  xBar = zeros( size (x) );
  if numel( z ) == 0, z = 0; end
  
  optIter = 0;
  relDiff = Inf;

  tauk = tau;
  sigmak = sigma;
  alpha0 = 0.95;
  eta = 0.95;
  c = 0.9;

  gamma = 2;
  rng(20240124)
  %%% data
  A = randn(20, 20);
  x = randn(20, 1);
  b = A*x;
  
  %%% functions
  f = @(x) norm(A*x - b, 2)^2;
  % proxf = @(y,t) proxL2Sq( y, 1, b, A );
  
  P = inv(A'*A + eye(size(A))); % compute inverse once and use thousands of times
  
  proxf = @(y, t) P * (y - b);
  
  Rg = @(x, gamma) 2 * proxg(x, gamma) - x;
  Rf = @(x, gamma) 2 * proxf(x, gamma) - x;

  alphak = alpha0;

  %%% Averaged operator iteration is of the form
  % z_k+1 = (1 - alpha) z_k + alpha S z_k

  phik = phi;
  pk = 0;
  zk = 0;
  while optIter < N
    optIter = optIter + 1;

    ak = phik(1);

    akp1 = proxf(ak - pk, gamma);

    Sphik = (theta / gamma)*applyA(2*akp1 - ak) + zk;
    Sphik = Sphik - 2*proxgConj(Sphik, theta/gamma);
    Sphik = Sphik * alpha;

    zkp1 = (1 - alpha)*zk + (theta / gamma)*applyA(akp1); %???
    zkp1 = zkp1 - (1 - alpha)*(theta/gamma)*applyA(ak) - alpha*Sphik;

    phikp1 = (1 - alpha)*phik + alpha*Sphik;

    pkp1 = phikp1(1) - akp1;


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
