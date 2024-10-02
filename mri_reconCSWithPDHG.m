
function out = mri_reconCSWithPDHG( kData, varargin )

  p = inputParser;
  p.addParameter( 'wavSplit', [], @isnumeric );
  p.parse( varargin{:} );
  wavSplit = p.Results.wavSplit;

  sImg = [ size( kData, 1 ) size( kData, 2 ) ];
  if numel( wavSplit ) == 0
    wavSplit = makeWavSplit( sImg );
  end
  mask = ( kData ~= 0 );

  unknownIndxs = find( mask == 0 );

  function out = A( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      x = zeros( sImg );
      x( unknownIndxs ) = in;
      img = fftshift2( ifft2( ifftshift2( x ) ) );
      out = wtDaubechies2( img, wavSplit );
    else
      WhIn = iwtDaubechies2( in, wavSplit );
      fftWhIn = fftshift2( ifft2h( ifftshift2( WhIn ) ) );
      out = fftWhIn( unknownIndxs );
    end
  end

  function out = innerProd( in1, in2 )   %#ok<DEFNU>
    out = real( dotP( in1, in2 ) );
  end

  doCheckAdjoint = false;
  if doCheckAdjoint
    [ checkA, errA ] = checkAdjoint( rand( numel( unknownIndxs ), 1 ), @A, ...
      'innerProd', @innerProd );   %#ok<UNRCH>
    if checkA == 0
      error([ 'A adjoint test failed with err: ', num2str( errA ) ]);
    else
      disp( 'A adjoint test passed' );
    end
  end

  function out = f( in )   %#ok<INUSD>
    out = 0;
  end

  function out = proxf( in, t )   %#ok<INUSD>
    out = in;
  end

  b = wtDaubechies2( fftshift2( ifft2( ifftshift2( kData ) ) ), wavSplit );
  function out = g( in )
    out = norm( in + b, 1 );
  end

  function out = proxgConj( in, sigma )
    out = in - sigma * proxL1Complex( in/sigma, 1/sigma, 1, -b );
  end

  normA = powerIteration( @A, rand( numel( unknownIndxs ), 1 ) );
  tau = 1d-1 / normA;
  [xStar,objValues,relDiffs] = pdhg( zeros( numel( unknownIndxs ), 1 ), @proxf, @proxgConj, tau, ...
    'A', @A, 'f', @f, 'g', @g, 'N', 10000, 'normA', normA, 'printEvery', 20, 'verbose', true, 'tol', [] );   %#ok<ASGLU>

  %[xStar,objValues] = pdhgWLS( zeros( numel( unknownIndxs ), 1 ), @proxf, @proxgConj, ...
  %  'N', 10000, 'A', @A, 'beta', 0.1, 'mu', 0.5, 'f', @f, 'g', @g, 'verbose', true );

figure;  semilogynice( objValues );

  kOut = kData;
  kOut( unknownIndxs ) = xStar;
  out = fftshift2( ifft2( ifftshift2( kOut ) ) );
end
