function testAdjointsAB2()

%% This is borked right now
n = 256;
sImg = [n n];
sFSR = [32 32];

idx = randi(prod(sImg), [17000, 1]);
id2 = sort(unique(idx));
nt = numel(id2);
wavSplit = makeWavSplit(sImg);

mask = zeros(n);
mask(id2) = 1;
lmask = mask == 1;


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
    out = kPF( id2 );
  end
  end

x = randn(n);
testx = randn([nt 1]);

[~, phaseimgx] = testA(x, sFSR);

phases = angle(phaseimgx);

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
Atx = ta(x, 'transp');

Bx = tb(testx);
Btx = tb(x, 'transp');

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

AAtx = ta(Atx);
BBtx = tb(Btx);

norm(AAtx + BBtx - 8*x / (n*n))
test1 = AAtx + BBtx;
test1 ./ x


end

