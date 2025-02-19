function output=ft2(in)
% Performs fftshift(fft2(ifftshift(input)))
% 2D forward FT
output = fftshift(fft2(ifftshift(in)));
