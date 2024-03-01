% load('kData_knee.mat', 'kData');
% load('brain.mat', 'd2');
% kData = d2;
load('ankle.mat')
kData = d1;

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


csDRRecons = cell(1,1,8);
csPDHGRecons = cell(1,1,8);
for coilIdx = 1:8
  
  coilData = fftSamples_wavACR(:,:,coilIdx);

  csDRRecons{1,1,coilIdx} = mri_reconCSWithDR( coilData, 'wavSplit', wavSplit );
  csPDHGRecons{1,1,coilIdx} = mri_reconCSWithPDHG( coilData, 'wavSplit', wavSplit );
end

csDRRecons = cell2mat( csDRRecons );  csDRRecon = mri_reconRoemer( csDRRecons );
figure; imshowscale(abs(csDRRecon)); title('cs with DR');

csPDHGRecons = cell2mat( csPDHGRecons );  csPDHGRecons = mri_reconRoemer( csPDHGRecons );
figure; imshowscale(abs(csPDHGRecons)); title('cs with PDHG');
