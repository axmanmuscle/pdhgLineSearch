%load('kData_knee.mat', 'kData');
load('brain.mat', 'd2');
kData = d2;
% load('ankle.mat')
% kData = d1;
rng(202402291);

kData = kData./max(abs(kData(:)));

gt_recon = mri_reconSSQ(kData);


sImg = size( kData, (1:2) );

vdSigFrac = 0.3;
sampleFraction = 0.3;

wavSplit = makeWavSplit( sImg );
[ fsr, sFSR ] = mri_makeFullySampledCenterRegion( sImg, wavSplit );

nSamples = round( sampleFraction * prod( sImg ) );
vdSig = round( vdSigFrac * sImg );

sampleMask = mri_makeSampleMask( sImg, nSamples, vdSig );
wavMaskACR = mri_makeSampleMask( sImg, sum(sampleMask(:)), vdSig, 'startMask', fsr>0 );

fftSamples_wavACR = bsxfun( @times, kData, wavMaskACR );

max_iters = 2500;

psnr_no_linesearch = zeros([max_iters, 1]);
psnr_linesearch = zeros([max_iters, 1]);

ssim_no_linesearch = zeros([max_iters, 1]);
ssim_linesearch = zeros([max_iters, 1]);

max_val = max(gt_recon, [], 'all');

nolsObjs = cell(8, 1);
lsObjs = cell(8,1);

for iters = max_iters:max_iters
    cspf_nolsRecons = cell(1,1,8);
    cspf_lsRecons = cell(1,1,8);
    for coilIdx = 1:8
    
      coilData = fftSamples_wavACR(:,:,coilIdx);
      fftSamples_wavACR_pf = coilData;
      fftSamples_wavACR_pf( ceil( sImg(1) / 2 ) + round( sFSR(1) / 2 ) : end, : ) = 0;
      fftSamples_wavACR_pf( fsr > 0 ) = coilData( fsr > 0 );
    
    
      [~,phaseImg] = mri_reconPartialFourier( fftSamples_wavACR_pf, sFSR );
      phases = angle( phaseImg );
      [nolsRecon, nolsObj] = psnr_plot_helper( fftSamples_wavACR_pf, sFSR, 'wavSplit', wavSplit, 'alg', 'primalDualDR_avgOp', 'N', iters );
      [lsRecon, lsObj] = psnr_plot_helper( fftSamples_wavACR_pf, sFSR, 'wavSplit', wavSplit, 'alg', 'primalDualDR_avgOp_wls', 'N', iters );
      cspf_nolsRecons{1,1,coilIdx} = nolsRecon;
      cspf_lsRecons{1,1,coilIdx} = lsRecon;

      nolsObjs{coilIdx} = nolsObj;
      lsObjs{coilIdx} = lsObj;
      
    end
    
    cspf_nolsRecons = cell2mat( cspf_nolsRecons ); noLinesearchRecon = mri_reconRoemer( cspf_nolsRecons );
    cspf_lsRecons = cell2mat( cspf_lsRecons ); linesearchRecon = mri_reconRoemer( cspf_lsRecons );

    psnr_nols = psnr(noLinesearchRecon, gt_recon, max_val);
      ssim_nols = ssim(noLinesearchRecon, gt_recon);
      psnr_no_linesearch(iters) = psnr_nols;
      ssim_no_linesearch(iters) = ssim_nols;

      psnr_ls = psnr(linesearchRecon, gt_recon, max_val);
      ssim_ls = ssim(linesearchRecon, gt_recon);
      psnr_linesearch(iters) = psnr_ls;
      ssim_linesearch(iters) = ssim_ls;
end

nolsObjsmat = cell2mat(nolsObjs);
lsObjsmat = cell2mat(lsObjs);
% figure; plot(nolsObjs); hold on; plot(lsObjs, 'bo'); title('objective values');
% figure; imshowscale(abs(linesearchRecon)); title('line search')

% absdiff = abs(linesearchRecon - cspfRecon);
% figure; imshowscale(absdiff);
