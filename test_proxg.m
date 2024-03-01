function test_proxg()

sImg = [4 4];
Amat = 4 * eye(4);
  
  function out = A(x, op)
    if nargin < 2 || strcmp( op, 'notransp' )
      out = Amat*x;
    else
      out = Amat'*x;
    end
  end

function out = A_full( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      img = fftshift2( ifft2( ifftshift2( in ) ) );
      out = wtDaubechies2( img );
    else
      WhIn = iwtDaubechies2( in );
      out = fftshift2( ifft2h( ifftshift2( WhIn ) ) );
    end
  end

  function out = g( in, lambda )
  out = lambda*norm( A_full( in ) , 1 );
  end

  function out = proxg( in, sigma, lambda )
  alpha = 16;
  out = in + alpha * A_full( proxL1Complex( A_full(in), lambda*sigma / alpha ) - A_full(in), 'transp') ;
  end

% function out = proxg( in, sigma, lambda )
%   alpha = 1/16;
%   out = in + alpha * A( proxL1Complex( A(in), lambda*sigma / alpha ) - A(in), 'transp') ;
%   end

  function out = h( in, lambda )
    out = lambda*norm(in, 1);
  end

  function out = proxh( in, t, lambda )
    out = softThresh(in, t*lambda);
  end

  function out = prox_def_g(x, v, t, lambda)
  out = g(x, lambda) + norm(x - v, 'fro')^2 / (2*t);
  end

  function out = prox_def_h(x, v, t, lambda)
  out = h(x, lambda) + norm(x - v, 'fro')^2 / (2*t);
  end


v = randn(sImg);
lambda = 6;
t = 6;

fun = @(x) prox_def_g(x, v, t, lambda);
% fun = @(x) prox_def_h(x, v, t, lambda);
options = optimoptions('fminunc', 'MaxFunctionEvaluations', 1e8, ...
                       'StepTolerance', 1e-10, 'OptimalityTolerance', 1e-12);%, ...
                       %'Display', 'iter-detailed');
% xtest = proxh(v, t, lambda);
xtest = proxg(v, t, lambda);
fxtest = fun(xtest);
s = 0;
for i = 1:50
  x0 = randn(sImg);
  xStar = fminunc(fun, x0, options);
  fxstar = fun(xStar);
  if fxstar < fxtest
    if abs(fxtest - fxstar) > 1e-9
      error('prox operator failed, optimization found smaller obj value')
    end
  end
  s = s + norm(xtest - xStar, 'fro');
end

fprintf('\nfinal diff is %f\n', s / 15);


end