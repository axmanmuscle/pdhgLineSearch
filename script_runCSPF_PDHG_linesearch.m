% load('kData_knee.mat', 'kData');
load('brain.mat', 'd2');
kData = d2;
% load('ankle.mat')
% kData = d1;
rng(202402291);

kData = kData./max(abs(kData(:)));

pfSamples = kData;

sImg = size( kData, (1:2) );

vdSigFrac = 0.3;
sampleFraction = 0.3;

wavSplit = makeWavSplit( sImg );
% wavSplit = zeros(2);  wavSplit(1,1) = 1;
[ fsr, sFSR ] = mri_makeFullySampledCenterRegion( sImg, wavSplit );

nSamples = round( sampleFraction * prod( sImg ) );
vdSig = round( vdSigFrac * sImg );

sampleMask = mri_makeSampleMask( sImg, nSamples, vdSig );
wavMaskACR = mri_makeSampleMask( sImg, sum(sampleMask(:)), vdSig, 'startMask', fsr>0 );

fftSamples_wavACR = bsxfun( @times, kData, wavMaskACR );


% csOut = mri_reconCSWithPDHG( fftSamples_wavACR, 'wavSplit', wavSplit );
csRecons = cell(1,1,8);
pfRecons = cell(1,1,8);
cspfRecons = cell(1,1,8);
linesearchRecons = cell(1,1,8);
parfor coilIdx = 1:8
% for coilIdx = 1:8
  pfData = pfSamples(:, :, coilIdx);
  coilData = fftSamples_wavACR(:,:,coilIdx);
  fftSamples_wavACR_pf = coilData;
%   fftSamples_wavACR_pf = pfData;
  % fftSamples_wavACR( ceil( ( sImg(1) + 1 ) / 2 ) + round( sFSR(1) / 2 ) : end, :, : ) = 0;
  fftSamples_wavACR_pf( ceil( sImg(1) / 2 ) + round( sFSR(1) / 2 ) : end, : ) = 0;
  fftSamples_wavACR_pf( fsr > 0 ) = pfData( fsr > 0 );

  [~,phaseImg] = mri_reconPartialFourier( fftSamples_wavACR_pf, sFSR );
  phases = angle( phaseImg );
%   pfRecons{1,1,coilIdx} = mri_reconPFHomodyne(fftSamples_wavACR_pf, sFSR, 'phases', phases);
  % pfRecons{1,1,coilIdx} = mri_reconPartialFourier( fftSamples_wavACR_pf, sFSR, 'phases', phases );
  % csRecons{1,1,coilIdx} = mri_reconCSWithPDHG( fftSamples_wavACR_pf, 'wavSplit', wavSplit );
  % cspfRecons{1,1,coilIdx} = mri_reconCSPFWithPDHG( fftSamples_wavACR_pf, sFSR, 'wavSplit', wavSplit  );
  %t1 = mri_reconCSPFHomodyne( fftSamples_wavACR_pf, sFSR, 'wavSplit', wavSplit  );
  % linesearchRecons{1,1,coilIdx} = mri_reconCSPFHomodyne_alex_pdhg( fftSamples_wavACR_pf, sFSR, 'wavSplit', wavSplit, 'alg', 'primalDualDR_avgOp_wls'  );
  % linesearchRecons{1,1,coilIdx} = mri_reconCSPFHomodyne_nick2( fftSamples_wavACR_pf, sFSR, 'wavSplit', wavSplit, 'alg', 'primalDualDR_avgOp_wls' );

  linesearchRecons{1,1,coilIdx} = mri_reconCSPFHomodyne_nick2( fftSamples_wavACR_pf, sFSR, 'wavSplit', wavSplit, 'alg', 'primalDualDR_avgOp' );

end
% csRecons = cell2mat( csRecons );  csRecon = mri_reconRoemer( csRecons );
% % pfRecons = cell2mat( pfRecons );  pfRecon = mri_reconRoemer( pfRecons );
% cspfRecons = cell2mat( cspfRecons );  cspfRecon = mri_reconRoemer( cspfRecons );
linesearchRecons = cell2mat( linesearchRecons ); linesearchRecon = mri_reconRoemer( linesearchRecons );

% figure; imshowscale(wavMaskACR);
% figure; imshowscale(abs(fftSamples_wavACR_pf) > 0);
% figure; imshowscale(abs(csRecon)); title('cs')
% figure; imshowscale(abs(pfRecon));  title('pf')
% figure; imshowscale(abs(cspfRecon)); title('cspf')
figure; imshowscale(abs(linesearchRecon)); title('line search')

% absdiff = abs(linesearchRecon - cspfRecon);
% figure; imshowscale(absdiff);
