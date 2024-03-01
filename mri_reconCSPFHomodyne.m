function recon =  mri_reconCSPFHomodyne(kData, sFSR, varargin)

p = inputParser;
p.addParameter( 'wavSplit', [], @isnumeric );
p.parse( varargin{:} );
wavSplit = p.Results.wavSplit;

sImg = [ size( kData, 1 ) size( kData, 2 ) ];
if numel( wavSplit ) == 0
  wavSplit = makeWavSplit( sImg );
end
mask = ( kData ~= 0 );
n0 = size(mask);

unknownMask = 1 - mask;
unknownMask( ceil(sImg(1)/2) + round( sFSR(1)/2 ) : end, : ) = 0;
unknownIndxs = find(unknownMask == 1);

[~,phaseImg] = mri_reconPFHomodyne( kData, sFSR );
phases = angle( phaseImg );

  function out = applyA(in, op)
    if nargin < 2 || strcmp( op, 'notransp' )
      x = zeros(n0);
      x(unknownIndxs) = in;
      Px = mri_reconPFHomodyne(x, sFSR, 'phases', phases);
      out = wtDaubechies2(Px, wavSplit);
    else
      Whin = iwtDaubechies2(in, wavSplit);
      Ptx = mri_reconPFHomodyne(Whin, sFSR, 'phases', phases, 'op', 'transp');
      out = Ptx(unknownIndxs);
    end
  end

  function out = applyB(in, op)

    Ny = size( in, 1 );
    ys = size2imgCoordinates( Ny );
    centerIndx = find( ys == 0, 1 );

    firstY = centerIndx - ceil( ( sFSR(1) - 1 ) / 2 );
    lastY = centerIndx + floor( ( sFSR(1) - 1 ) / 2 );
    m = 2 / ( lastY - firstY );
    ramp = m * ys + 1.0;
    ramp( 1 : firstY ) = 0;
    ramp( lastY : end ) = 2;
    ramp = 2 - ramp;
    ramp2 = ramp.*ramp;
    ramp3 = bsxfun(@times, mask, ramp2);
    alpha = 1/8;
    nramp = size(ramp);
    bramp = (1/alpha) *ones(nramp) - ramp3;
    sbramp = sqrt(bramp);

    if nargin < 2 || strcmp( op, 'notransp' )
      Py = mri_reconPFHomodyne(in, sFSR, 'phases', phases, 'ramp', sbramp);
      out = wtDaubechies2(Py, wavSplit);
    else
      WhIn = iwtDaubechies2(in, wavSplit);
      out = mri_reconPFHomodyne(WhIn, sFSR, 'phases', phases, 'ramp', sbramp, 'op', 'transp');
    end
  end

  function out = proxftilde( in, t )   %#ok<INUSD>
    out = zeros(size(in));
    out(1:n) = in(1:n);
  end

  function out = proxL1(in, sigma, b)
    out = softThresh(in + b, sigma) - b;
  end

  PsiPb = wtDaubechies2( mri_reconPFHomodyne( kData, sFSR, 'phases', phases ), wavSplit ) ;
  function out = proxgconj( in, sigma )
    out = in - sigma * proxL1( in/sigma, 1/sigma, PsiPb );
  end

  function out = g(in)
    out = norm(in + PsiPb, 1);
  end

  function out = gtilde(in)
    x = in(1:n);
    y = in(n+1:end);
    y = reshape(y, n0);
    out = g( applyA(x) + applyB(y) );
  end


  function out = ftilde(in)
%     y = in(n+1:end);
%     if norm(y, 'inf') > 0
%       out = Inf;
%     else
%       out = 0;
%     end
    out = 0;
  end

  function out = proxgtilde(in, tau)
    %%% slide 12.8 Vandenberghe's notes 236C
    sigma = theta / tau;
    x = in(1:n);
    y = in(n+1:end);
    y = reshape(y, n0);
    t = proxgconj( sigma * ( applyA(x) + applyB(y) ), sigma);
    Att = applyA( t, 'transp' );
    Btt = applyB( t, 'transp' );
    out = in - tau * [Att; Btt(:)];
  end

  function out = proxgtildeconj( in, lambda )
    out = in - lambda*proxgtilde(in/lambda, 1/lambda);
  end

  function out = proxf( in, t )
    out = in;
  end

  function out = f(in)
    out = 0;
  end

testx = randn(size(kData));
ATk = applyA(testx, 'transp');
AATk = applyA(ATk);

BTk = applyB(testx, 'transp');
BBTk = applyB(BTk);

test1 = AATk + BBTk;
test1 ./ testx


img0 = fftshift2( ifft2( ifftshift2( kData ) ) );
x0 = wtDaubechies2( img0, wavSplit );
value = findValueBelowFraction( x0, 0.9 );
x0( x0 < value ) = 0;
img0 = iwtDaubechies2( x0, wavSplit );
fftImg0 = fftshift2( fft2( ifftshift2( img0 ) ) );
k0 = fftImg0( unknownIndxs );
n = numel(k0);
z0 = zeros(n0);

theta = 8/prod(sImg);

x0 = [k0; z0(:)];
maxIter = 10;

normA = powerIteration(@applyA, rand( numel( unknownIndxs ), 1 ) );
tau = 0.95 / normA;
sigma_pdhg = 0.95 / normA;

% [xStar,objValues,relDiffs] = pdhg( k0, @proxf, @proxgconj, tau, ...
%     'sigma', sigma_pdhg, 'A', @applyA, 'f', @f, 'g', @g, 'N', 1000, 'normA', normA, 'printEvery', 50, ...
%     'verbose', true, 'tol', [] );

% xStar = douglasRachford(x0, @proxftilde, @proxgtilde, 1, 'verbose', true);
% xStar = primal_dual_dr_2(x0, @proxftilde, @proxgtildeconj);
xStar = primal_dual_dr_aoi_wls(x0, @proxftilde, @proxgtildeconj, @ftilde, @gtilde, maxIter);

out = kData;
out(unknownIndxs) = xStar(1:n);

recon = mri_reconPFHomodyne(out, sFSR, 'phases', phases);

end