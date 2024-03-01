
function out = mri_reconCSWithDR( kData, varargin )

p = inputParser;
p.addParameter( 'wavSplit', [], @isnumeric );
p.parse( varargin{:} );
wavSplit = p.Results.wavSplit;

sImg = [ size( kData, 1 ) size( kData, 2 ) ];
if numel( wavSplit ) == 0
  wavSplit = makeWavSplit( sImg );
end
mask = ( kData ~= 0 );

knownIndxs = find( mask == 1 );

  function out = A( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      img = fftshift2( ifft2( ifftshift2( in ) ) );
      out = wtDaubechies2( img, wavSplit );
    else
      WhIn = iwtDaubechies2( in, wavSplit );
      out = fftshift2( ifft2h( ifftshift2( WhIn ) ) );
    end
  end

  function out = innerProd( in1, in2 )   %#ok<DEFNU>
    out = real( dotP( in1, in2 ) );
  end

doCheckAdjoint = true;
if doCheckAdjoint
  [ checkA, errA ] = checkAdjoint( rand( size( kData ) ), @A, ...
    'innerProd', @innerProd );   %#ok<UNRCH>
  if checkA == 0
    error([ 'A adjoint test failed with err: ', num2str( errA ) ]);
  else
    disp( 'A adjoint test passed' );
  end
end

b = kData(knownIndxs);

  function out = f( in )   %#ok<INUSD>
    if norm(in(knownIndxs) - b) == 0
      out = 0;
    else
      out = Inf;
    end
  end

  function out = proxf( in, t )   %#ok<INUSD>
    out = in;
    out( knownIndxs ) = b;
  end

%   b = wtDaubechies2( fftshift2( ifft2( ifftshift2( kData ) ) ), wavSplit );
  function out = g( in, lambda )
    out = lambda*norm( A( in ) , 1 );
  end

  function out = proxg( in, sigma  )
    alpha = prod(sImg);
    out = in + alpha * A( proxL1Complex( A(in), sigma / alpha ) - A(in), 'transp') ;
  end

lambda = prod(sImg);
pg = @(x, t) proxg(x, t*lambda);

x0 = zeros(size(kData));

[xStar,objValues,alphas] = dr_wLS( x0, @proxf, pg, 10);   %#ok<ASGLU>

out = fftshift2( ifft2( ifftshift2( xStar ) ) );

end
