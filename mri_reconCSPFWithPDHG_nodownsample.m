
function out = mri_reconCSPFWithPDHG_nodownsample( kData, sFSR, varargin )

  p = inputParser;
  p.addParameter( 'wavSplit', [], @isnumeric );
  p.parse( varargin{:} );
  wavSplit = p.Results.wavSplit;

  sImg = [ size( kData, 1 ) size( kData, 2 ) ];
  if numel( wavSplit ) == 0
    wavSplit = makeWavSplit( sImg );
  end
  mask = ( kData ~= 0 );
  knownIndxs = find( mask == 1);

  unknownMask = 1 - mask;
  unknownMask( ceil(sImg(1)/2) + round( sFSR(1)/2 ) : end, : ) = 0;
  unknownIndxs = find( unknownMask == 1 );

  [~,phaseImg] = mri_reconPartialFourier( kData, sFSR );
  phases = angle( phaseImg );

  function out = A( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
%       x = zeros( sImg );
%       x( unknownIndxs ) = in;
      img = mri_reconPartialFourier( in, sFSR, 'phases', phases );
      out = wtDaubechies2( img, wavSplit );
    else
      WhIn = iwtDaubechies2( in, wavSplit );
      out = mri_reconPartialFourier( WhIn, sFSR, 'phases', phases, 'op', 'transp' );
%       out = kPF( unknownIndxs );
    end
  end


  function out = innerProd( in1, in2 )   %#ok<DEFNU>
    out = real( dotP( in1, in2 ) );
  end

  doCheckAdjoint = true;
  if doCheckAdjoint
    [ checkA, errA ] = checkAdjoint( rand( size(kData)), @A, ...
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
    out = zeros(size(in));
    out(knownIndxs) = kData(knownIndxs);
  end

%   b = wtDaubechies2( mri_reconPartialFourier( kData, sFSR, 'phases', phases ), wavSplit );
  function out = g( in )
    out = norm( in, 1 );
  end

  function out = proxgConj( in, sigma )
    out = in - sigma * proxL1Complex( in/sigma, 1/sigma, 1);
  end

  img0 = fftshift2( ifft2( ifftshift2( kData ) ) );
  x0 = wtDaubechies2( img0, wavSplit );
  value = findValueBelowFraction( x0, 0.9 );
  x0( x0 < value ) = 0;
  img0 = iwtDaubechies2( x0, wavSplit );
  fftImg0 = fftshift2( fft2( ifftshift2( img0 ) ) );
%   k0 = fftImg0( unknownIndxs );
k0 = fftImg0;

  normA = powerIteration( @A, rand( size(kData) ) );
  tau = 0.95d-1 / normA;
  sigma = 0.95d1 / normA;

  [xStar,objValues,relDiffs] = pdhg( zeros(size(kData)), @proxf, @proxgConj, tau, ...
    'sigma', sigma, 'A', @A, 'f', @f, 'g', @g, 'N', 1000, 'normA', normA, 'printEvery', 50, ...
    'verbose', true, 'tol', [] );   %#ok<ASGLU>

  %[xStar,objValues] = pdhgWLS( zeros( numel( unknownIndxs ), 1 ), @proxf, @proxgConj, ...
  %  'N', 10000, 'A', @A, 'beta', 2, 'mu', 0.5, 'f', @f, 'g', @g, 'printEvery', 25, 'verbose', true );

  kOut = kData;
  kOut( unknownIndxs ) = xStar;
  out = mri_reconPartialFourier( kOut, sFSR, 'phases', phases );
end
