
function [x,objValues] = avgOpIter_wLS( x0, S_in, varargin )
% Implements an averaged opterator iteration with line search.
% See "Line Search for Averaged Operator Iteration" by Gisellson
% et al. (2016)
%
% x = avgOpIter( x0, S [, 'alpha', alpha, 'maxIter', maxIter ] )
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
p.addParameter( 'N', 100, @ispositive );
p.addParameter( 'objFunction', [] );
p.addParameter( 'verbose', false );
p.parse( varargin{:} );
N = p.Results.N;
objFunction = p.Results.objFunction;
verbose = p.Results.verbose;

if nargout > 2
    if numel( objFunction ) == 0
        error( 'Must specify an objective function to return objective values' );
    end
    objValues = zeros( N, 1 );
end

%%% parameters
alpha_bar = 0.5; % alpha_bar
eps = 0.03; % eps for (1 - eps) || rbar_k || in linesearch
tol = 1e-7; % tolerance for exit criterion
alpha_change = 1/1.4; % factor for change in alpha during linesearch

if isnumeric(S_in)
    S = @(x) S_in*x;
else
    S = @(x) S_in(x);
end

x = x0;
for optIter = 1 : N
    rk = S(x) - x;
    x_bar = x + alpha_bar*rk;
    rk_bar = S(x_bar) - x_bar;

    alpha_k = 50;
    subiter = 0;
    while true
        subiter = subiter + 1;
        % fprintf('subiter %d\n', subiter);
        xkp1 = x + alpha_k*rk;
        rkp1 = S(x) - x;
        if norm(rkp1) < (1-eps)*norm(rk_bar)
            x = xkp1;
            break
        end
        alpha_k = alpha_k*alpha_change;
        if alpha_k < alpha_bar
            alpha_k = alpha_bar;
            xkp1 = x + alpha_k*rk;
            x = xkp1;
            break
        end
    end

    if nargout > 1 || ( numel(objFunction) > 0 && verbose == true )
        objValue = objFunction( x );
    end

    if nargout > 1
        objValues( optIter ) = objValue;
    end

    if verbose == true
        outStr = [ 'avgOpIter: Completed ', indx2str(optIter,N), ' of ', num2str(N) ];
        if numel( objFunction ) > 0
            outStr = [ outStr, '  objective: ', num2str( objValue ) ];   %#ok<AGROW>
        end
        disp( outStr );
    end

end

end
