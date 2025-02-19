%	parameters for the 3 disk phantom
%	x,y center, radius, 'amplitude' (e.g. attenuation coefficient)
circ = [0 0 75 0.1; 30 30 18 0.2; -53 0 9 0.3];
nobj = size(circ,1);

%	image parameters
nx = 192; ny = 192;
dx = 1;		                 % 1 mm / pixel

%	geometry parameters
nr = 192;	dr = 1;		% # of radial samples, ray spacing

na = round(pi * nr/2);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 1 - determine number of angular samples
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = dr * [-nr/2:nr/2-1]';

function ang = ang_pos(na)
    ang = [0:(na-1)]'/na * pi;	% angular sample positions
    fprintf('Number of angles = %g\n', na)
end

ang = ang_pos(na);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 2 - compute sinogram for disk phantom (NO CHANGES NEEDED)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sinogram1 = GenerateSinogram(nr,na,nobj,circ,ang,r)
    rr = r(:,ones(1,na));
    sinogram1 = zeros(nr, na);
    for ii=1:nobj
        cx = circ(ii,1);		% center of object in x
        cy = circ(ii,2);		% center of object in y
        rad = circ(ii,3);		% radius of object
        amp = circ(ii,4);		% amplitude of object (in atten
    
        % correct amplitude for overlying objects
        if ii > 1, amp = amp - circ(1,4);, end
    
        % center location of object for each projection
        tau = cx * cos(ang) + cy * sin(ang);
        tau = tau(:,ones(1,nr))';
    
        % find all locations where "rr" is within "rad" of "tau"
        t = find( (rr-tau).^2 <= rad.^2 );
    
        % update the sinogram with length of segment (a bit of geometry...)
        sinogram1(t) = sinogram1(t)+amp*2*sqrt(rad^2-(rr(t)-tau(t)).^2);
    
    end
end

sinogram1 = GenerateSinogram(nr,na,nobj,circ,ang,r);


% Output Image of Singram
figure(1)
imagesc(r,ang,sinogram1'); colormap('gray')
title('Sinogram of Disk Phantom')
xlabel('Detector Position (mm)')
ylabel('Angular Position (radians)')

xlim([-93 94])
ylim([0.02 3.08])

sinogram = sinogram1;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 3 - Implement plain backprojection
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bpimage = SimpleBPImage(nx,ny,na,sinogram,nr,ang, verbose)
    bpimage = zeros(nx,ny);
    for ia = 1:na
        % fprintf('Angle %g of %g\n', ia, na)
    
        % first backproject at theta = 0 (copy the contents of the projection
        % across the image)
        tmpim = repmat(sinogram(:,ia), 1, nr);
    
        % now rotate the projection
        rotated_proj = imrot3(tmpim, ang(ia), 'bilinear');
    
        % Add projection to the image
        bpimage = bpimage + rotated_proj;

        if verbose
            figure(ia)
            x = (-nx/2:nx/2) * dx % x-coordinate (mm)
            y = (-ny/2:ny/2) * dx % y-coordinate (mm)
            imagesc(x,y,bpimage'); colormap('gray'); axis('image');axis('xy');
            title('Simple Backprojection Image')
            xlabel('X - Spatial Coordinate (mm)')
            ylabel('Y - Spatial Coordinate (mm)')
        end
    
    end
end
bpimage = SimpleBPImage(nx,ny,na,sinogram,nr,ang,0);

% Display Image
figure(2)
x = (-nx/2:nx/2) * dx % x-coordinate (mm)
y = (-ny/2:ny/2) * dx % y-coordinate (mm)
imagesc(x,y,bpimage'); colormap('gray'); axis('image');axis('xy');
title('Simple Backprojection Image')
xlabel('X - Spatial Coordinate (mm)')
ylabel('Y - Spatial Coordinate (mm)')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 4 - Filter Projections
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sinogramfilt = FilterProjection(sinogrampad)
    % Apply Fourier Transform
    sinogramft = ft(sinogrampad);
    
    % Create the Ramp Filter
    freqs = linspace(-0.5, 0.5, size(sinogramft, 1))';
    RampFilter = abs(freqs);  
    
    % Apply the Ramp Filter
    FilteredFFT = sinogramft .* RampFilter;  
    
    % Inverse FFT to get filtered projections
    sinogramfilt = real(ift(FilteredFFT));
end


% zero pad sinogram (if using Fourier methods)
sinogrampad = padarray(sinogram, [nr/2, 0], 'both');

sinogramfilt = FilterProjection(sinogrampad);

% Trim padding
sinogramfilt_trimmed = sinogramfilt(nr/2+1:nr+nr/2, :);


%
% Plot Filtered Sinogram
%

figure(3)
plot(r, sinogram(:,1)./max(sinogram(:,1)), '-',...
    r, sinogramfilt_trimmed(:,1)./max(sinogramfilt_trimmed(:,1)),':');
title('Projection at theta = 0')
xlabel('Detector Position (mm)')
ylabel('Normalized Intensity')
legend('before filtering','after filtering', location ='best')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 5 - Backproject the filtered sinogram
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imagefbp = SimpleBPImage(nx,ny,na,sinogramfilt_trimmed,nr,ang,0);

%
% Display Reconstructed Image with Negatives Set to Zero
%

figure(4)
x = (-nx/2:nx/2-1) * dx;  % X-coordinates in mm
y = (-ny/2:ny/2-1) * dx;  % Y-coordinates in mm
imagesc(x, y, max(imagefbp', 0)); colormap('gray'); axis('image'); axis('xy');
title('FBP Reconstruction Image')
xlabel('Spatial Coordinate - X (mm)')
ylabel('Spatial Coordinate - Y (mm)')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 6 - Generate Fourier Interpolation Image
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to generate Fourier Interpolation Image
function output=FI_image(in, nx, ny, nr, dr, ang)
    fov = dr*nr;
    kx = 1/fov.*[-nx/2:nx/2-1];
    ky = 1/fov.*[-ny/2:ny/2-1];
    kr = 1/fov.*[-nr/2:nr/2-1];
    [kaa,krr] = meshgrid(-ang,kr);
    [kxx,kyy] = meshgrid(kx,ky);
    kxxin = -krr.*sin(kaa);
    kyyin = krr.*cos(kaa);
    
    sinogramft = ft(in);  % 1D of FT of projections
    fprintf('kxxin %g x %g,kyyin %g x %g,sinogramft %g x %g,kxx %g x %g,kyy %g x %g', ...
        size(kxxin),size(kyyin),size(sinogramft),size(kxx),size(kyy))
    fdata = griddata(kxxin,kyyin,sinogramft,kxx,kyy,'linear');
    t = find(isnan(fdata));
    fdata(t) = zeros(size(t));
    imagefi = ift2(fdata);
    
    output = imagefi;
end

image_fi = FI_image(sinogram, nx, ny, nr, dr, ang)

%
% Display Reconstructed Image
%
figure(5)
x = (-nx/2:nx/2-1) * dx;  % X-coordinates in mm
y = (-ny/2:ny/2-1) * dx;  % Y-coordinates in mm
imagesc(x, y, abs(image_fi')); colormap('gray'); axis('image'); axis('xy')
title('FI Reconstruction Image')
xlabel('Spatial Coordinate - X (mm)')
ylabel('Spatial Coordinate - Y (mm)')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 7 - Generate Fourier Interpolation Image using zeropadding
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Zero-pad the sinogram in the R dimension to 384 samples.
nr_zp = 384
nx_zp = 384
ny_zp = 384
padding_amount = (nr_zp - nr) / 2
sinogram_zp = padarray(sinogram,[padding_amount, 0], 'both')

% FI image with more samples
imagefizp = FI_image(sinogram_zp, nx_zp,ny_zp,nr_zp,dr,ang)

% Crop to only focuson the middle 192x192 image
imagefizp_cropped = real(imagefizp(nx_zp/4:3*nx_zp/4 - 1, ny_zp/4:3*ny_zp/4 -1));

figure(6)
x = (-nx/2:nx/2-1) * dx;  % X-coordinates in mm
y = (-ny/2:ny/2-1) * dx;  % Y-coordinates in mm
imagesc(x, y, abs(imagefizp_cropped')); colormap('gray'); axis('image'); axis('xy');
title('FI Reconstruction Image with ZP')
xlabel('Spatial Coordinate - X (mm)')
ylabel('Spatial Coordinate - Y (mm)')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 8 - Plot profiles through reconstructed images
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


y = dx * ([1:ny]'-(ny+1)/2);
center_row = 133;%floor(ny/2) + 1;
linefbp = imagefbp(center_row, :);
linefi  = abs(image_fi(center_row, :));  % take absolute value (magnitude)
linefizp = abs(imagefizp_cropped(center_row, :));

figure(7)
plot(x, linefbp./max(linefbp), '-', x, linefi./max(linefi), ':', x, linefizp/max(linefizp), '--');
xlabel('Radial Position (mm)')
ylabel('Normalized Projection Value')
legend('FBP image', 'FI image', 'Oversampled FI image', Location='best')
title('Profile through Reconstruction (y = 0)')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 9 - Do again for subsampled image
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Sub Sampled Sinogram %%%%%%%
sub_na = floor(na/4)
sub_ang = ang_pos(sub_na)
sinosubsamp = GenerateSinogram(nr,sub_na,nobj,circ,sub_ang,r)

% Output Image of Singram
figure(8)
imagesc(r,sub_ang,sinosubsamp'); colormap('gray')
title('Subsampled Sinogram of Disk Phantom')
xlabel('Detector Position (mm)')
ylabel('Angular Position (radians)')


%%%%%%%%%%%%%% FBP %%%%%%%%%%%%%%%%%
% Padding Sinogram
sinogrampad = padarray(sinosubsamp, [nr/2, 0], 'both');
% Filter Projection
sinogramfilt = FilterProjection(sinogrampad);
% Trim padding
sinogramfilt_trimmed = sinogramfilt(nr/2+1:nr+nr/2, :);
% Generating image
imagefbp = SimpleBPImage(nx,ny,sub_na,sinogramfilt_trimmed,nr,sub_ang,0);

% Display Reconstructed Image with Negatives Set to Zero
figure(9)
x = (-nx/2:nx/2-1) * dx;  % X-coordinates in mm
y = (-ny/2:ny/2-1) * dx;  % Y-coordinates in mm
imagesc(x, y, max(imagefbp', 0)); colormap('gray'); axis('image'); axis('xy');
title('FBP Reconstruction Sub-sampled Image')
xlabel('Spatial Coordinate - X (mm)')
ylabel('Spatial Coordinate - Y (mm)')

%%%%%%%%%%%%%%% FI %%%%%%%%%%%%%%%%%
image_fi = FI_image(sinosubsamp, nx, ny, nr, dr, sub_ang)

% Display Reconstructed Image
figure(10)
x = (-nx/2:nx/2-1) * dx;  % X-coordinates in mm
y = (-ny/2:ny/2-1) * dx;  % Y-coordinates in mm
imagesc(x, y, abs(image_fi')); colormap('gray'); axis('image'); axis('xy')
title('FI Reconstruction Sub-sampled Image')
xlabel('Spatial Coordinate - X (mm)')
ylabel('Spatial Coordinate - Y (mm)')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 10 - Do again for lmited view angles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ang = ang_pos(na)
quarter_ang = ang(1:sub_na)
sinosubsamp = GenerateSinogram(nr,sub_na,nobj,circ,quarter_ang,r)

% Output Image of Singram
figure(11)
imagesc(r,sub_ang,sinosubsamp'); colormap('gray')
title('Subsampled Sinogram of Disk Phantom (Limited View Angles)')
xlabel('Detector Position (mm)')
ylabel('Angular Position (radians)')


%%%%%%%%%%%%%% FBP %%%%%%%%%%%%%%%%%
% Padding Sinogram
sinogrampad = padarray(sinosubsamp, [nr/2, 0], 'both');
% Filter Projection
sinogramfilt = FilterProjection(sinogrampad);
% Trim padding
sinogramfilt_trimmed = sinogramfilt(nr/2+1:nr+nr/2, :);
% Generating image
imagefbp = SimpleBPImage(nx,ny,sub_na,sinogramfilt_trimmed,nr,quarter_ang,0);

% Display Reconstructed Image with Negatives Set to Zero
figure(12)
x = (-nx/2:nx/2-1) * dx;  % X-coordinates in mm
y = (-ny/2:ny/2-1) * dx;  % Y-coordinates in mm
imagesc(x, y, max(imagefbp', 0)); colormap('gray'); axis('image'); axis('xy');
title('FBP Reconstruction Sub-sampled Image (Limited View Angles)')
xlabel('Spatial Coordinate - X (mm)')
ylabel('Spatial Coordinate - Y (mm)')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Question 11 - load and reconstruct mystery object
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load mys25; % gives ang and sinogram
na = length(ang);
figure(13)
imagesc(r,ang,sinogram'); colormap('gray')
title('Sinogram of Mystery Object')
xlabel('Detector Position (mm)')
ylabel('Angular Position (radians)')

%%%%%%%%%%%%%% FBP %%%%%%%%%%%%%%%%%
% Padding Sinogram
sinogrampad = padarray(sinogram, [nr/2, 0], 'both');
% Filter Projection
sinogramfilt = FilterProjection(sinogrampad);
% Trim padding
sinogramfilt_trimmed = sinogramfilt(nr/2+1:nr+nr/2, :);
% Generating image
imagefbp = SimpleBPImage(nx,ny,na,sinogramfilt_trimmed,nr,ang,0);

% Display Reconstructed Image with Negatives Set to Zero
figure(14)
x = (-nx/2:nx/2-1) * dx;  % X-coordinates in mm
y = (-ny/2:ny/2-1) * dx;  % Y-coordinates in mm
imagesc(x, y, max(imagefbp', 0)); colormap('gray'); axis('image'); axis('xy');
title('FBP Reconstruction Image of Mystery Object')
xlabel('Spatial Coordinate - X (mm)')
ylabel('Spatial Coordinate - Y (mm)')