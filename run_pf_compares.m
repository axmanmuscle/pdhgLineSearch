load('kData_knee.mat', 'kData');
% load('brain.mat', 'd2');
% kData = d2;
% load('ankle.mat')
% kData = d1;

kData = kData./max(abs(kData(:)));

sImg = size( kData, (1:2) );

gt_recon = mri_reconSSQ(kData);

vdSigFrac = 0.3;
sampleFraction = 0.3;

wavSplit = makeWavSplit( sImg );
[ fsr, sFSR ] = mri_makeFullySampledCenterRegion( sImg, wavSplit );

nSamples = round( sampleFraction * prod( sImg ) );
vdSig = round( vdSigFrac * sImg );

sampleMask = mri_makeSampleMask( sImg, nSamples, vdSig );
wavMaskACR = mri_makeSampleMask( sImg, sum(sampleMask(:)), vdSig, 'startMask', fsr>0 );

fftSamples_wavACR = bsxfun( @times, kData, wavMaskACR );


nu = 5/8;
pf_fs_idx = ceil(sImg(1) * nu);
sFSR_fs = sImg(1)*(nu - 1/2) * 2;

csRecons = cell(1,1,8);
pfRecons = cell(1,1,8);
cspfRecons = cell(1,1,8);

pf_fs_Recons = cell(1,1,8);
% parfor coilIdx = 1:8
for coilIdx = 1:8
  
  coilData = fftSamples_wavACR(:,:,coilIdx);
  fftSamples_wavACR_pf = coilData;
  fftSamples_wavACR_pf( ceil( sImg(1) / 2 ) + round( sFSR(1) / 2 ) : end, : ) = 0;
  fftSamples_wavACR_pf( fsr > 0 ) = coilData( fsr > 0 );

  coilData_fs = kData(:, :, coilIdx);
  coilData_fs( pf_fs_idx : end, : ) = 0;

  [~,phaseImg] = mri_reconPartialFourier( fftSamples_wavACR_pf, sFSR );
  phases = angle( phaseImg );
  pfRecons{1,1,coilIdx} = mri_reconPartialFourier( fftSamples_wavACR_pf, sFSR, 'phases', phases );
  csRecons{1,1,coilIdx} = mri_reconCSWithPDHG( fftSamples_wavACR_pf, 'wavSplit', wavSplit );
  cspfRecons{1,1,coilIdx} = mri_reconCSPFWithPDHG( fftSamples_wavACR_pf, sFSR, 'wavSplit', wavSplit  );

  [~, phaseImg_fs] = mri_reconPartialFourier(coilData_fs, sFSR_fs);
  phases_fs = angle(phaseImg_fs);
  pf_fs_Recons{1,1,coilIdx} = mri_reconPartialFourier( coilData_fs, sFSR, 'phases', phases_fs );
end
pfRecons = cell2mat( pfRecons );  pfRecon = mri_reconRoemer( pfRecons );
csRecons = cell2mat( csRecons );  csRecon = mri_reconRoemer( csRecons );
cspfRecons = cell2mat( cspfRecons ); cspfRecon = mri_reconRoemer( cspfRecons ); 
pf_fs_Recons = cell2mat( pf_fs_Recons );  pf_fs_recon = mri_reconRoemer( pf_fs_Recons );

figure; imshowscale(abs(pfRecon)); title('pf')
figure; imshowscale(abs(csRecon)); title('cs')
figure; imshowscale(abs(cspfRecon)); title('cspf')
figure; imshowscale(abs(pf_fs_recon)); title('pf_fs')

max_gt = max(gt_recon, [], 'all');
[psnr_pf, snr_pf] = psnr(abs(pfRecon), gt_recon, max_gt);
ssim_pf = ssim(abs(pfRecon, gt_recon));

[psnr_cs, snr_cs] = psnr(abs(csRecon), gt_recon, max_gt);
ssim_cs = ssim(abs(csRecon, gt_recon));

[psnr_cspf, snr_cspf] = psnr(abs(cspfRecon), gt_recon, max_gt);
ssim_cspf = ssim(abs(cspfRecon, gt_recon));

% absdiff = abs(cspfRecon - csRecon);
% figure; imshowscale(absdiff);
