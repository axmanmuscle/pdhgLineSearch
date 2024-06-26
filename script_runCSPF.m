% load('kData_knee.mat', 'kData');
load('brain.mat', 'd2');
kData = d2;
% load('ankle.mat')
% kData = d1;

kData = kData./max(abs(kData(:)));

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
%parfor coilIdx = 1:8
for coilIdx = 1:8
  
  coilData = fftSamples_wavACR(:,:,coilIdx);
  fftSamples_wavACR_pf = coilData;
  fftSamples_wavACR_pf( ceil( sImg(1) / 2 ) + round( sFSR(1) / 2 ) : end, : ) = 0;
  fftSamples_wavACR_pf( fsr > 0 ) = coilData( fsr > 0 );

  [~,phaseImg] = mri_reconPartialFourier( fftSamples_wavACR_pf, sFSR );
  phases = angle( phaseImg );
%   pfRecons{1,1,coilIdx} = mri_reconPartialFourier( fftSamples_wavACR_pf, sFSR, 'phases', phases );
%   csRecons{1,1,coilIdx} = mri_reconCSWithPDHG( fftSamples_wavACR_pf, 'wavSplit', wavSplit );
%   cspfRecons{1,1,coilIdx} = mri_reconCSPFWithPDHG( fftSamples_wavACR_pf, sFSR, 'wavSplit', wavSplit  );
%   t1 = mri_reconCSPFWithPDHG( fftSamples_wavACR_pf, sFSR, 'wavSplit', wavSplit  );
  t2 = mri_reconCSPFWithPDHG_nodownsample( fftSamples_wavACR_pf, sFSR, 'wavSplit', wavSplit  );
end
csRecons = cell2mat( csRecons );  csRecon = mri_reconRoemer( csRecons );
pfRecons = cell2mat( pfRecons );  pfRecon = mri_reconRoemer( pfRecons );
cspfRecons = cell2mat( cspfRecons );  cspfRecon = mri_reconRoemer( cspfRecons );

% figure; imshowscale(wavMaskACR);
% figure; imshowscale(abs(fftSamples_wavACR_pf) > 0);
figure; imshowscale(abs(csRecon)); 
figure; imshowscale(abs(pfRecon)); 
figure; imshowscale(abs(cspfRecon)); 

absdiff = abs(cspfRecon - csRecon);
figure; imshowscale(absdiff);
