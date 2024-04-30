%% partial fourier recon

load('kData_knee.mat', 'kData');

kData_orig = kData;

if_kdata = mri_reconIFFT(kData, 'multislice', true);
image_ssq = mri_reconRoemer(if_kdata);
kData = ifftshift2(fft2(fftshift2(image_ssq)));

n = size(kData, 1);
mid = n/2 + 1;

kData_test = kData;
kData_test(mid+1:end, :) = 0;

kData_test(mid+1:end, 1) = conj(flip(kData_test(2:mid-1, 1)));

flip2 = kData_test(2:mid-1, 2:end);
kData_test(mid+1:end, 2:end) = conj(flip(flip(flip2, 1), 2));

norm(kData - kData_test)
