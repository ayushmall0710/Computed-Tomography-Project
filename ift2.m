function output=ift2(in)
% Performs fftshift(ifft2(ifftshift( input)))
% 2D inverse FT
output = fftshift(ifft2(ifftshift(in)));
