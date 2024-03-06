
function [x,objValues,alphas] = avgOpIter_wLS_fast( x0, S_in, varargin )
  % Implements an averaged opterator iteration with line search.
  % See "Line Search for Averaged Operator Iteration" by Gisellson et al. (2016)
  %
  % x = avgOpIter( x0, S [, 'alpha_bar', alpha_bar, 'maxIter', maxIter ] )
  %
  % Inputs:
  % x0 - the initial guess
  % S - either a matrix or a a function handle that is the non-expansive operator
  %
  % Written by Nicholas - Copyright 2024
  %
  % https://github.com/ndwork/dworkLib.git
  %
  % This software is offered under the GNU General Public License 3.0.  It
  % is offered without any warranty expressed or implied, including the
  % implied warranties of merchantability or fitness for a particular
  % purpose.
  
  p = inputParser;
  p.addParameter( 'alpha_bar', 0.5, @(x) x>=0 && x<=1 );
  p.addParameter( 'N', 100, @ispositive );
  p.addParameter( 'objFunction', [] );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'verbose', false );
  p.parse( varargin{:} );
  alpha_bar = p.Results.alpha_bar;
  N = p.Results.N;
  objFunction = p.Results.objFunction;
  printEvery = p.Results.printEvery;
  verbose = p.Results.verbose;
  
  if nargout > 1
      if numel( objFunction ) == 0
          error( 'Must specify an objective function to return objective values' );
      end
      objValues = zeros( N, 1 );
  end
  
  if nargout > 2
    alphas = zeros( N, 1 );
  end
  
  %%% parameters
  eps = 0.01; % eps for (1 - eps) || rbar_k || in linesearch
  alpha0 = 2;
  %alpha_change = 1/1.4; % factor for change in alpha during linesearch
  alpha_change = 0.95; % factor for change in alpha during linesearch

  if isnumeric( S_in )
      S = @(x) S_in*x;
  else
      S = @(x) S_in(x);
  end
  
  x = x0;
  rk = S(x) - x;
  
  for optIter = 1 : N
    x_bar = x + alpha_bar * rk;
    rk_bar = S(x_bar) - x_bar;
    normRkBar = sqrt( real( dotP( rk_bar, rk_bar ) ) );

    alpha_k = alpha0;
    subiter = 0;
    while true
      subiter = subiter + 1;
      xkp1 = x + alpha_k * rk;
      rkp1 = S( xkp1 ) - xkp1;
      normRkp1 = sqrt( real( dotP( rkp1, rkp1 ) ) );

      if normRkp1 < (1-eps) * normRkBar
        x = xkp1;
        rk = rkp1;
        break
      end

      alpha_k = alpha_k * alpha_change;

      if alpha_k <= alpha_bar
        alpha_k = alpha_bar;
        x = x_bar;
        rk = rk_bar;
        break
      end
    end

    if nargout > 1 || ( numel(objFunction) > 0 && verbose == true )
      objValue = objFunction( x );
    end

    if nargout > 1
      objValues( optIter ) = objValue;
    end
    if nargout > 2
      alphas( optIter ) = alpha_k;
    end

    if verbose == true && mod( optIter, printEvery ) == 0
      outStr = [ 'avgOpIter: Completed ', indx2str(optIter,N), ' of ', num2str(N), ...
        '   alpha: ', num2str(alpha_k) ];
      if numel( objFunction ) > 0
          outStr = [ outStr, '  objective: ', num2str( objValue ) ];   %#ok<AGROW>
      end
      disp( outStr );
    end

  end

end
