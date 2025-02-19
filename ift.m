function output=ift(in)
% Performs ifftshift(ifft(ifftshift( input)))
% 1D inverse FT
output = fftshift(ifft(ifftshift(in)));
