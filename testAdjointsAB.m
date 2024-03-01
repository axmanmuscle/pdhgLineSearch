function testAdjointsAB()

n = 256;
sFSR = [32 32];

x = randn(n);
y = randn(n);

% x = randn(n) + 1i*randn(n);
% y = randn(n) + 1i*randn(n);

[~, phaseimgx] = testA(x, sFSR);
[~, phaseimgy] = testA(y, sFSR);

phasex = angle(phaseimgx);
phasey = angle(phaseimgy);

dotprod1 = @(x, y) real( y(:)' * x(:));

Ax = testA(x, sFSR, 'phases', phasex);
Ay = testA(y, sFSR, 'phases', phasey);

Atx = testA(x, sFSR, 'phases', phasey, 'op', 'transp');
Aty = testA(y, sFSR, 'phases', phasex, 'op', 'transp');

Bx = testB(x, sFSR, 'phases', phasex);
By = testB(y, sFSR, 'phases', phasey);

Btx = testB(x, sFSR, 'phases', phasey, 'op', 'transp');
Bty = testB(y, sFSR, 'phases', phasex, 'op', 'transp');

%%% <Ax, y> = <x, Aty>
t1 = dotprod1(Ax, y);
t2 = dotprod1(x, Aty);
fprintf('Ax vs Aty test: %f\n',t1 - t2);

%%% <x, Ay> = <Atx, y>
t3 = dotprod1(x, Ay);
t4 = dotprod1(Atx, y);
fprintf('Ay vs Atx test: %f\n',t3 - t4);

%%% <Ax, x> = <x, Atx>
%%% have to recompute Atx with x's phase
[Atx, Atxphase] = testA(x, sFSR, 'phases', phasex, 'op', 'transp');
t1 = dotprod1(Ax, x);
t2 = dotprod1(x, Atx);
fprintf('Ax vs Atx test: %f\n',t1 - t2);

%%% <Ay, y> = <y, Aty>
%%% have to recompute Aty with y's phase
Aty = testA(y, sFSR, 'phases', phasey, 'op', 'transp');
t1 = dotprod1(Ay, y);
t2 = dotprod1(y, Aty);
fprintf('Ay vs Aty test: %f\n',t1 - t2);

%%% <Bx, y> = <x, Bty>
t1 = dotprod1(Bx, y);
t2 = dotprod1(x, Bty);
fprintf('Bx vs Bty test: %f\n',t1 - t2);

%%% <x, By> = <Btx, y>
t3 = dotprod1(x, By);
t4 = dotprod1(Btx, y);
fprintf('By vs Btx test: %f\n',t3 - t4);

%%% <Bx, x> = <x, Btx>
%%% have to recompute Btx
[Btx, Btxphase] = testB(x, sFSR, 'phases', phasex, 'op', 'transp');
t1 = dotprod1(Bx, x);
t2 = dotprod1(x, Btx);
fprintf('Bx vs Btx test: %f\n',t1 - t2);

[AAtx, AAtxphase] = testA(Atx, sFSR, 'phases', phasex);
[BBtx, BBtxphase] = testB(Btx, sFSR, 'phases', phasex);

norm(AAtx + BBtx - 8*x / (n*n))
test1 = AAtx + BBtx;
test1 ./ x


end

function ax = Aop(x, op)
  global zphase
  sFSR = [8, 8];
  ax = testA(x, sFSR, 'phases', zphase, 'op', op);
end
