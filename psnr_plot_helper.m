function [recon, objValues] =  psnr_plot_helper(kData, sFSR, varargin )
  % recon =  mri_reconCSPFHomodyne(kData, sFSR [, 'wavSplit', wavSplit ] )

  p = inputParser;
  p.addParameter( 'alg', 'douglasRachford' );
  p.addParameter( 'doChecks', false );
  p.addParameter( 'gamma', 1d-2, @ispositive );
  p.addParameter( 'printEvery', 1, @ispositive );
  p.addParameter( 'wavSplit', [], @isnumeric );
  p.addParameter( 'N', 100, @isnumeric);
  p.parse( varargin{:} );
  alg = p.Results.alg;
  doChecks = p.Results.doChecks;
  gamma = p.Results.gamma;
  printEvery = p.Results.printEvery;
  wavSplit = p.Results.wavSplit;

  iters = p.Results.N;
  
  sImg = [ size( kData, 1 ) size( kData, 2 ) ];
  if numel( wavSplit ) == 0
    wavSplit = makeWavSplit( sImg );
  end
  mask = ( kData ~= 0 );
  sMask = size( mask );

  unknownMask = 1 - mask;
  unknownMask( ceil( ( sImg(1) + 1 ) / 2 ) + round( sFSR(1)/2 ) : end, : ) = 0;
  unknownIndxs = find( unknownMask == 1 );
  nUnknown = numel( unknownIndxs );

  [~,phaseImg] = mri_reconPFHomodyne( kData, sFSR );
  phases = angle( phaseImg );

  function out = applyA(in, op)
    if nargin < 2 || strcmp( op, 'notransp' )
      x = zeros(sMask);
      x( unknownIndxs ) = in;
      Px = mri_reconPFHomodyne( x, sFSR, 'phases', phases);
      out = wtDaubechies2( Px, wavSplit );
    else
      Whin = iwtDaubechies2(in, wavSplit);
      Ptx = mri_reconPFHomodyne( Whin, sFSR, 'phases', phases, 'op', 'transp' );
      out = Ptx( unknownIndxs );
    end
  end

  Ny = sMask(1);
  function out = applyB(in, op)
    ys = size2imgCoordinates( Ny );
    centerIndx = find( ys == 0, 1 );

    firstY = centerIndx - ceil( ( sFSR(1) - 1 ) / 2 );
    lastY = centerIndx + floor( ( sFSR(1) - 1 ) / 2 );
    m = 2 / ( lastY - firstY );
    ramp = m * ys + 1.0;
    ramp( 1 : firstY ) = 0;
    ramp( lastY : end ) = 2;
    ramp = 2 - ramp;
    ramp2 = ramp .* ramp;
    ramp3 = bsxfun(@times, unknownMask, ramp2);
    rampB = sqrt( numel( mask ) * ones( sMask ) - ramp3 );
    % This scaling inside the square root makes A A^T + B B^T = I

    in = reshape( in, sMask );
    if nargin < 2 || strcmp( op, 'notransp' )
      Py = mri_reconPFHomodyne(in, sFSR, 'phases', phases, 'ramp', rampB);
      out = wtDaubechies2(Py, wavSplit);
    else
      WhIn = iwtDaubechies2(in, wavSplit);
      out = mri_reconPFHomodyne(WhIn, sFSR, 'phases', phases, 'ramp', rampB, 'op', 'transp');
    end
  end

  nMask = prod( sMask );
  function out = applyAB( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      out = applyA( in(1:nUnknown) ) + applyB( in(nUnknown+1:end) );
    else
      out = zeros( nUnknown + nMask, 1 );
      out(1:nUnknown) = applyA( in, 'transp' );
      out(nUnknown+1:end) = applyB( in, 'transp' );
    end
  end

  if doChecks == true
    innerProd = @(x,y) real( dotP( x, y ) );

    x0 = rand( numel( unknownIndxs ), 1 );
    [chkA, errChkA] = checkAdjoint( x0, @applyA, 'innerProd', innerProd );
    if chkA == false
      error(['Check for adjoint of A failed with error ', num2str(errChkA) ]);
    else
      disp( 'Check for adjoint of A passed.' );
    end

    y0 = rand( sMask );
    [chkB, errChkB] = checkAdjoint( y0, @applyB, 'innerProd', innerProd );
    if chkB == false
      error(['Check for adjoint of B failed with error ', num2str(errChkB) ]);
    else
      disp( 'Check for adjoint of B passed.' );
    end

    % Check to make sure A A^T + B B^T is a scaled identity
    AATpBBTy = applyA( applyA( y0, 'transp' ) ) + applyB( applyB( y0, 'transp' ) );
    testAB = AATpBBTy ./ y0;
    errAB = std( testAB(:) );
    if errAB > 1d-12
      error([ 'AB check failed with error ', num2str(errAB) ]);
    else
      disp( 'AB check passed' );
    end
  end

  function out = f( in )   %#ok<INUSD>
    out = 0;
  end

  function out = proxf( in, t )   %#ok<INUSD>
    out = in;
  end

  function out = f_tilde( in )
    if max( abs( in( nUnknown + 1 : end ) ) ) > 0
      out = Inf;
    else
      out = 0;
    end
  end

  function out = proxf_tilde( in, t )   %#ok<INUSD>
    out = zeros( size( in ) );
    out( 1 : nUnknown ) = in( 1 : nUnknown );
  end

  PsiPMcTb = wtDaubechies2( mri_reconPFHomodyne( kData, sFSR, 'phases', phases ), wavSplit ) ;
  function out = g( in )
    out = norm( in(:) + PsiPMcTb(:), 1 );
  end

  function out = proxg( in, t )
    out = softThresh( in + PsiPMcTb, t ) - PsiPMcTb;
  end

  function out = proxgConj( in, t )
    out = proxConj( @proxg, in, t );
  end

  function out = g_tilde( in )
    x = in( 1 : nUnknown );
    y = reshape( in( nUnknown + 1 : end ), sMask );
    AxBy = applyA( x ) + applyB( y );
    out = norm( AxBy(:) + PsiPMcTb(:), 1 );
  end

  function out = proxg_tilde( in, t )
    out = proxCompositionAffine( @proxg, in, @applyAB, 0, 1, t );
  end

  function out = proxg_tildeConj( in, t )
    out = proxConj( @proxg_tilde, in, t );
  end

  function out = Rf( in, gamma )
    out = 2 * proxf_tilde( in, gamma ) - in;
  end

  function out = Rg( in, gamma )
    out = 2 * proxg_tilde( in, gamma ) - in;
  end

  function out = S_DR( in )
    out = Rg( Rf( in, gamma ), gamma );
  end

  function out = Rg_tildeConj( in, gamma )
    out = 2 * proxg_tildeConj( in, gamma ) - in;
  end

  function out = S_pdDR( in )
    out = -gamma * Rg_tildeConj( Rf( in, gamma ) / gamma , 1/gamma );
  end

  k0 = zeros( nUnknown, 1 );
  z0 = zeros( sMask );
  x0 = [ k0; z0(:) ];

  %normA = [];
  %if numel( gamma ) == 0
  %  normA = powerIteration( @applyA, rand( size( k0 ) ) );
  %  gamma = 1 / normA;
  %end

  switch alg

    case 'douglasRachford'
      [xStar,objValuesDR] = douglasRachford( x0, @proxf_tilde, @proxg_tilde, gamma, 'N', 200, ...
        'f', @f_tilde, 'g', @g_tilde, 'verbose', true, 'printEvery', printEvery );   %#ok<ASGLU>

    case 'douglasRachford_avgOp'
      objF = @(x) f( proxf_tilde(x) ) + g( proxf_tilde(x) );
      [xStar] = avgOpIter( x0, @S_DR, 'alpha', 0.5, 'N', 1000, ...
        'objFunction', objF, 'verbose', true, 'printEvery', printEvery );

    case 'pdhg'
      tau = gamma;
      if numel( normA ) == 0, normA = powerIteration( @applyA, rand( size( k0 ) ) ); end
      [kStar,objValues] = pdhg( k0, @proxf, @proxgConj, tau, 'A', @applyA, 'normA', normA, 'sigma', gamma, ...
        'N', 1000, 'f', @f, 'g', @g, 'verbose', true, 'printEvery', printEvery );
      xStar = zeros( size( x0 ) );
      xStar( 1 : nUnknown ) = kStar;

    case 'primalDualDR'
      [xStar, objValues] = primalDualDR( x0, @proxf_tilde, @proxg_tildeConj, gamma, 'N', 1000, ...
        'f', @f_tilde, 'g', @g_tilde, 'verbose', true);%, 'printEvery', printEvery );   %#ok<ASGLU>

    case 'primalDualDR_avgOp'
      pddrAvgOpObjF = @(x) g_tilde( proxf_tilde( x ) );
      [xStar,objValues] = avgOpIter( x0, @S_pdDR, 'alpha', 0.5, 'N', iters, ...
        'objFunction', pddrAvgOpObjF, 'verbose', true, 'printEvery', 100 );

    case 'primalDualDR_avgOp_wls'
      objF = @(x) f_tilde( proxf_tilde(x) ) + g_tilde( proxf_tilde(x) );
      [xStar,objValues,alphas] = avgOpIter_wLS( x0, @S_pdDR, 'N', 100, ...
        'objFunction', objF, 'verbose', true, 'printEvery', 100, 'doLineSearchTest', true );   %#ok<ASGLU>

    otherwise
      error( 'Unrecognized algorithm' );
  end

  kOut = kData;
  kOut( unknownIndxs ) = xStar( 1 : nUnknown );
  recon = mri_reconPFHomodyne( kOut, sFSR, 'phases', phases);
end
