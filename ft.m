function output=ft(in)
% Performs fftshift(fft(ifftshift(input)))
% 1-D forward FT
output = fftshift(fft(ifftshift(in)));
