function testAdjointsAB2()

rng(20240229);
n = 64;
sImg = [n n];
sFSR = [16 16];

idx = randi(prod(sImg), [17000, 1]);
id2 = sort(unique(idx));
unknownIndxs = id2;
% load('ukinxs.mat', 'unknownIndxs')
% id2 = unknownIndxs;
nt = numel(id2);
wavSplit = makeWavSplit(sImg);

mask = zeros(n);
mask(id2) = 1;
lmask = mask == 1;

n0 = n;


  function out = ta(in, op)
  if nargin < 2 || strcmp( op, 'notransp' )
    tx = zeros( sImg );
    tx( id2 ) = in;
    img = testA( tx, sFSR, 'phases', phases );
    out = wtDaubechies2( img, wavSplit );
  else
    WhIn = iwtDaubechies2( in, wavSplit );
    kPF = testA( WhIn, sFSR, 'phases', phases, 'op', 'transp' );
    out = kPF( id2 );
  end
  end

  function out = tb( in, op )
  if nargin < 2 || strcmp( op, 'notransp' )
    img = testB_2( in, sFSR, mask, 'phases', phases );
    out = wtDaubechies2( img, wavSplit );
  else
    WhIn = iwtDaubechies2( in, wavSplit );
    kPF = testB_2( WhIn, sFSR, mask, 'phases', phases, 'op', 'transp' );
    out = kPF;
  end
  end

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

x1 = randn(n);
testx = randn([nt 1]);

[~, phaseimgx] = testA(x1, sFSR);

phases1 = angle(phaseimgx);
phases = phases1;

dotprod1 = @(x, y) real( y(:)' * x(:));


% x = randn([1024 1]);
% y = randn([256*256 1]);

  function out = L(in)
    out = zeros([n*n 1]);
    out(id2) = in;
  end

Lt = @(in) in(id2);

% dotprod1(A(x), y)
% dotprod1(x, A(y, 'transp'))

Ax = ta(testx);
Atx1 = ta(x1, 'transp');

Bx = tb(x1);
Btx1 = tb(x1, 'transp');

% %%% <Ax, x> = <x, Atx>
% %%% have to recompute Atx with x's phase
% t1 = dotprod1(Ax, x);
% t2 = dotprod1(testx, Atx);
% fprintf('Ax vs Atx test: %f\n',t1 - t2);
% 
% %%% <Bx, x> = <x, Btx>
% %%% have to recompute Btx
% t1 = dotprod1(Bx, x);
% t2 = dotprod1(x, Btx);
% fprintf('Bx vs Btx test: %f\n',t1 - t2);

% [AAtx, AAtxphase] = testA(Atx, sFSR, 'phases', phasex);
% [BBtx, BBtxphase] = testB(Btx, sFSR, 'phases', phasex);

AAtx1 = ta(Atx1);
BBtx1 = tb(Btx1);

norm(AAtx1 + BBtx1 - 8*x1 / (n*n))
test1 = AAtx1 + BBtx1;
test1 ./ x1

Atx2 = applyA(x1, 'transp');
Btx2 = applyB(x1, 'transp');

AAtx2 = ta(Atx2);
BBtx2 = tb(Btx2);

AAtx3 = applyA(Atx2);
BBtx3 = applyB(Btx2);

test2 = AAtx2 + BBtx2;
test3 = AAtx3 + BBtx3;

end

