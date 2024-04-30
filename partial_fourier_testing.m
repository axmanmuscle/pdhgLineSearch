%% partial fourier idk


n = 256;
rng(20240416);
b = randn(n);
% fb = ifftshift2(fft2(fftshift2(b)));
fb = fft2(b);
mid = n/2 + 1;

fb_test = fb;
fb_test(mid+1:end, :) = 0;

fb_test(mid+1:end, 1) = conj(flip(fb_test(2:mid-1, 1)));

flip2 = fb_test(2:mid-1, 2:end);
fb_test(mid+1:end, 2:end) = conj(flip(flip(flip2, 1), 2));

norm(fb - fb_test)
