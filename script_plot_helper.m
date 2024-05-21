% load('kData_knee.mat', 'kData');
% load('brain.mat', 'd2');
% kData = d2;
load('ankle.mat')
kData = d1;
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

iters = 10000;

max_val = max(gt_recon, [], 'all');

nolsObjs = cell(8, 1);
lsObjs = cell(8,1);
pdhgObjs = cell(8,1);

cspf_nolsRecons = cell(1,1,8);
cspf_lsRecons = cell(1,1,8);
cspf_pdhgRecons = cell(1,1,8);
csRecons = cell(1,1,8);
gamma = 2^-10;
pdhg_gam = 2^10;
parfor coilIdx = 1:8

    coilData = fftSamples_wavACR(:,:,coilIdx);
    fftSamples_wavACR_pf = coilData;
    fftSamples_wavACR_pf( fsr > 0 ) = coilData( fsr > 0 );
    fftSamples_wavACR_pf( ceil( sImg(1) / 2 ) + round( sFSR(1) / 2 ) : end, : ) = 0;
    
    [~,phaseImg] = mri_reconPartialFourier( fftSamples_wavACR_pf, sFSR );
    phases = angle( phaseImg );
    csRecons{1,1,coilIdx} = mri_reconCSWithPDHG( fftSamples_wavACR_pf, 'wavSplit', wavSplit );
    % [nolsRecon, nolsObj] = psnr_plot_helper( fftSamples_wavACR_pf, sFSR, 'wavSplit', wavSplit, 'alg', 'primalDualDR_avgOp', 'N', iters, 'gamma', gamma  );
    [lsRecon, lsObj] = psnr_plot_helper( fftSamples_wavACR_pf, sFSR, 'wavSplit', wavSplit, 'alg', 'primalDualDR_avgOp_wls', 'N', iters, 'gamma', gamma  );
    [pdhgRecon, pdhgObj] = psnr_plot_helper( fftSamples_wavACR_pf, sFSR, 'wavSplit', wavSplit, 'alg', 'pdhg', 'N', 2000, 'gamma', pdhg_gam  );
    % cspf_nolsRecons{1,1,coilIdx} = nolsRecon;
    cspf_lsRecons{1,1,coilIdx} = lsRecon;
    cspf_pdhgRecons{1,1,coilIdx} = pdhgRecon;

    % nolsObjs{coilIdx} = nolsObj;
    lsObjs{coilIdx} = lsObj;
    pdhgObjs{coilIdx} = pdhgObj;

end

csRecons = cell2mat(csRecons); csRecon = mri_reconRoemer(csRecons);
% cspf_nolsRecons = cell2mat( cspf_nolsRecons ); noLinesearchRecon = mri_reconRoemer( cspf_nolsRecons );
cspf_lsRecons = cell2mat( cspf_lsRecons ); linesearchRecon = mri_reconRoemer( cspf_lsRecons );
cspf_pdhgRecons = cell2mat( cspf_pdhgRecons ); pdhgRecon = mri_reconRoemer( cspf_pdhgRecons );

% psnr_nols = psnr(noLinesearchRecon, gt_recon, max_val);
% ssim_nols = ssim(noLinesearchRecon, gt_recon);

psnr_ls = psnr(linesearchRecon, gt_recon, max_val);
ssim_ls = ssim(linesearchRecon, gt_recon);

psnr_pdhg = psnr(pdhgRecon, gt_recon, max_val);
ssim_pdhg = ssim(pdhgRecon, gt_recon);

% nolsObjsmat = cell2mat(nolsObjs);
lsObjsmat = cell2mat(lsObjs);
pdhgObjsmat = cell2mat(pdhgObjs);
% figure; plot(nolsObjs); hold on; plot(lsObjs, 'bo'); title('objective values');
% figure; imshowscale(abs(linesearchRecon)); title('line search')

% absdiff = abs(linesearchRecon - cspfRecon);
% figure; imshowscale(absdiff);
