
function out = mri_reconCSPFWithPDHG_WLS( kData, sFSR, varargin )

p = inputParser;
p.addParameter( 'wavSplit', [], @isnumeric );
p.parse( varargin{:} );
wavSplit = p.Results.wavSplit;

sImg = [ size( kData, 1 ) size( kData, 2 ) ];
if numel( wavSplit ) == 0
  wavSplit = makeWavSplit( sImg );
end
mask = ( kData ~= 0 );
n0 = size(mask, 1);

unknownMask = 1 - mask;
unknownMask( ceil(sImg(1)/2) + round( sFSR(1)/2 ) : end, : ) = 0;
unknownIndxs = find( unknownMask == 1 );

[~,phaseImg] = mri_reconPartialFourier( kData, sFSR );
phases = angle( phaseImg );

  function out = A( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      x = zeros( sImg );
      x( unknownIndxs ) = in;
      img = testA( x, sFSR, 'phases', phases );
      out = wtDaubechies2( img, wavSplit );
    else
      WhIn = iwtDaubechies2( in, wavSplit );
      kPF = testA( WhIn, sFSR, 'phases', phases, 'op', 'transp' );
      out = kPF( unknownIndxs );
    end
  end

% function out = B( in, op )
%     if nargin < 2 || strcmp( op, 'notransp' )
%       x = zeros( sImg );
%       x( unknownIndxs ) = in;
%       img = testB( x, sFSR, 'phases', phases );
%       out = wtDaubechies2( img, wavSplit );
%     else
%       WhIn = iwtDaubechies2( in, wavSplit );
%       kPF = testB( WhIn, sFSR, 'phases', phases, 'op', 'transp' );
%       out = kPF( unknownIndxs );
%     end
%   end

  function out = B( in, op )
    if nargin < 2 || strcmp( op, 'notransp' )
      img = testB_2( in, sFSR, unknownMask, 'phases', phases );
      out = wtDaubechies2( img, wavSplit );
    else
      WhIn = iwtDaubechies2( in, wavSplit );
      kPF = testB_2( WhIn, sFSR, unknownMask, 'phases', phases, 'op', 'transp' );
      out = kPF;
    end
  end

  function out = innerProd( in1, in2 )   %#ok<DEFNU>
    out = real( dotP( in1, in2 ) );
  end

doCheckAdjoint = true;
if doCheckAdjoint
  [ checkA, errA ] = checkAdjoint( rand( numel( unknownIndxs ), 1 ), @A, ...
    'innerProd', @innerProd );   %#ok<UNRCH>
  if checkA == 0
    error([ 'A adjoint test failed with err: ', num2str( errA ) ]);
  else
    disp( 'A adjoint test passed' );
  end

  [ checkB, errB ] = checkAdjoint( rand( size(kData)), @B, ...
    'innerProd', @innerProd );   %#ok<UNRCH>
  if checkB == 0
    error([ 'B adjoint test failed with err: ', num2str( errB ) ]);
  else
    disp( 'B adjoint test passed' );
  end
end

  function out = gtilde(in)
    x = in(1:n);
    y = in(n+1:end);
    y = reshape(y, [n0 n0]);
    out = g(A(x) + B(y));
  end

  function out = ftilde(in)
    y = in(n+1:end);
    if norm(y, 'inf') > 0
      out = Inf;
    else
      out = 0;
    end
  end

  function out = proxftilde( in, t )   %#ok<INUSD>
    out = zeros(size(in));
    out(1:n) = in(1:n);
  end


b = wtDaubechies2( abs(mri_reconPartialFourier( kData, sFSR, 'phases', phases )), wavSplit );
  function out = g( in )
    out = norm( in + b, 1 );
  end

  function out = proxL1(in, sigma, b)
    out = softThresh(in - b, sigma) + b;
  end

  function out = proxgconj( in, sigma )
    out = in - sigma * proxL1( in/sigma, 1/sigma, -b );
  end

  function out = proxGtilde(in, tau)
    sigma = alpha / tau;
    tmp1 = in(n+1:end);
    tmp1 = reshape(tmp1, [n0 n0]);
    tmp1 = proxgconj(sigma*(applyA(in(1:n)) + applyB(tmp1)), sigma);
    tmp2 = applyAt(tmp1);
    tmp3 = applyBt(tmp1);
    out = in - tau * [tmp2; tmp3(:)];
  end

  function out = proxGTildeConj(in, lambda)
    out = in - lambda * proxGtilde(in/lambda, 1/lambda);
  end

%     function out = proxGtildestar(in, tau)
%         tmp1 = in(n+1:end);
%         tmp1 = reshape(tmp1, [n0 n0]);
%         sig = tau/alpha;
%         tmp1 = proxgconj(sig*(applyA(in(1:n)) + applyB(tmp1)), sig);
%         tmp2 = applyAt(tmp1);
%         tmp3 = applyBt(tmp1);
%         out = tau * [tmp2; tmp3(:)];
%     end

% load('testmat', 'x', 'Atx1', 'AAtx1', 'Btx1', 'BBtx1', 'phases1');
% testx = x;
% clear x;
% [~, phaseImgx] = testA(testx, sFSR);
% phases = angle(phaseImgx);
testx = randn(size(kData));

Atx = A(testx, 'transp');
AAtx = A(Atx);

Btx = B(testx, 'transp');
BBtx = B(Btx);

t1 = AAtx + BBtx;
t1 ./ testx;

img0 = fftshift2( ifft2( ifftshift2( kData ) ) );
x0 = wtDaubechies2( img0, wavSplit );
value = findValueBelowFraction( x0, 0.9 );
x0( x0 < value ) = 0;
img0 = iwtDaubechies2( x0, wavSplit );
fftImg0 = fftshift2( fft2( ifftshift2( img0 ) ) );
k0 = fftImg0( unknownIndxs );
n = numel(k0);
z0 = zeros(n0);

x0 = [k0; z0(:)];

normA = powerIteration( @A, rand( numel( unknownIndxs ), 1 ) );
%   tau = 0.1;
alpha = 8/prod(sImg);
%   sigma = tau/theta;
tau = 0.95d-1 / normA;
% sigma = 0.95d1 / normA;

applyA = @(x) A(x, 'notransp');
applyAt = @(x) A(x, 'transp');
applyB = @(x) B(x, 'notransp');
applyBt = @(x) B(x, 'transp');

maxIter = 100;

%   pgtilde = @(x, t) x - tau * [A'; B'] * pgstar(sigma*(applyA(x(1:n)) + applyB(x(n+1:end))), sigma);
%   [xStar, iters] = primal_dual_dr_3(k0, @proxf, @proxgconj, @A, @B, tau);
%   [xStar2, iters2] = primal_dual_dr_aoi([k0;zeros(size(k0))], @proxf, @proxGtilde);
%   [xStar3, iters3] = primal_dual_dr_aoi([k0;zeros(size(k0))], @proxf, pgtildestar);
%   [xStar4, iters4, alphas] = dr_wLS([k0;zeros(size(k0))], @proxf, @proxGtilde, maxIter);

[xStar5, iters5, alphas5, objVals] = primal_dual_dr_aoi_wls(x0, @proxftilde, @proxGTildeConj, @ftilde, @gtilde);

% [xStar,objValues,relDiffs] = pdhg( k0, @proxf, @proxgconj, tau, ...
%     'sigma', sigma, 'A', @A, 'f', @f, 'g', @g, 'N', 1000, 'normA', normA, 'printEvery', 50, ...
%     'verbose', true, 'tol', [] );   %#ok<ASGLU>

% kOut = kData;
% kOut( unknownIndxs ) = xStar;
% out1 = mri_reconPartialFourier( kOut, sFSR, 'phases', phases );
% %
%   kOut = kData;
%   kOut( unknownIndxs ) = xStar2(1:n);
%   out2 = mri_reconPartialFourier( kOut, sFSR, 'phases', phases );
%
%   kOut = kData;
%   kOut( unknownIndxs ) = xStar3(1:n);
%   out3 = mri_reconPartialFourier( kOut, sFSR, 'phases', phases );
%
%   kOut = kData;
%   kOut( unknownIndxs ) = xStar4(1:n);
%   out4 = mri_reconPartialFourier( kOut, sFSR, 'phases', phases );
plot(objVals);
kOut = kData;
kOut( unknownIndxs ) = xStar5(1:n);
out5 = mri_reconPartialFourier( kOut, sFSR, 'phases', phases );

% figure; imshowscale(abs(out1)); title('xStar')
%   figure; imshowscale(abs(out2)); title('xStar2')
%   figure; imshowscale(abs(out3)); title('xStar3')
%   figure; imshowscale(abs(out4)); title('xStar4')
%   figure; imshowscale(abs(out5)); title('xStar5')

out = abs(out5);
figure; imshowscale(out)
end
