% load('brain.mat', 'd2');
% kData_coils = d2;
load('kData_knee.mat', 'kData');
kData_coils = kData;

kData_coils = kData_coils./max(abs(kData_coils(:)));

if_kdata = mri_reconIFFT(kData_coils, 'multislice', true);
image_ssq = mri_reconRoemer(if_kdata);
kData = ifftshift2(fft2(fftshift2(image_ssq)));

sImg = size( kData, (1:2) );

nu = 5/8;

coilIdx = 1;
coilData = kData_coils(:, :, coilIdx);

fftSamples_wavACR_pf = coilData;
fftSamples_wavACR_pf( ceil( sImg(1) * nu )+1 : end, : ) = 0;

sFSR = sImg(1)*(nu - 1/2) * 2;
% sFSR = 32;

[img,phaseImg] = mri_reconPartialFourier( fftSamples_wavACR_pf, sFSR );
phases = angle( phaseImg );

