clear; close all; clc;

%% INPUTS

projectionName = input('Enter the name of the projected data: ','s');
data = load(projectionName);
projections = data.pt;
M = max(max(max(data.rows_teta)));
N = max(max(max(data.cols_teta)));

[numberOfProjections, numberOfSamples] = size(projections);

stepSize = (180 - 0) / numberOfProjections;
thetaValues = linspace(0, (180 - stepSize), numberOfProjections);

xCoor = floor(-M/2):ceil(M/2);
yCoor = floor(-N/2):ceil(N/2);

%% FILTER

projectionsFiltered = fRamlak(projections);

%% BACKPROJECTION

image = zeros(M,N);
imageFiltered = zeros(M,N);
distances = zeros(numberOfSamples, (length(xCoor)-1));
rows = zeros(numberOfSamples, (length(xCoor)-1));
cols = zeros(numberOfSamples, (length(yCoor)-1));
pt = zeros(numberOfSamples);
ptFiltered = zeros(numberOfSamples);

for teta = 1:numberOfProjections
    distances = data.distances_teta(:,:,teta);
    rows = data.rows_teta(:,:,teta);
    cols = data.cols_teta(:,:,teta);
    pt = projections(teta, :);
    ptFiltered = projectionsFiltered(teta, :);

    for t = 1:numberOfSamples
        for i = 1:(length(xCoor)-1)
            if((0 < rows(t,i)) && (0 < cols(t,i)))
                image(rows(t,i), cols(t,i)) = image(rows(t,i), cols(t,i)) + pt(t) * distances(t,i);
                imageFiltered(rows(t,i), cols(t,i)) = imageFiltered(rows(t,i), cols(t,i)) + ptFiltered(t) * distances(t,i);
            end
        end
    end
end

%% VISUALIZATION

load lena
normalized = lena / max(lena(:));
imageNormalized = image / max(image(:));
imageFilteredNormalized = imageFiltered / max(imageFiltered(:));

%[R,xp] = radon(normalized,thetaValues);
I = iradon(projections', thetaValues','linear','none');
I = I / max(I(:));
I2 = iradon(projections', thetaValues');
I2 = I2 / max(I2(:));

% figure;
% imshow(normalized);
% title('Original Image')
figure;
imshow(I);
title('iradon unfiltered')
figure;
imshow(I2);
title('iradon with filter')
figure;
imshow(imageNormalized);
title('BP Image unfiltered')
figure;
imshow(imageFilteredNormalized);
title('BP Image with filter')

%% RAMLAK
function sinogram = fRamlak(sinogram)
% Filter image for before backprojection in frequency domain
% using ramlak triangular filter
% 
% sinogram: image converted to radon space
% make a Ram-Lak filter. it's just abs(f).

sinogram = sinogram';
N1 = size(sinogram,1);
freqs=linspace(-1, 1, N1).';
myFilter = abs(freqs);
%myFilter = repmat(myFilter, [1 180]);

% FT domain filtering
ft_R = fftshift(fft(sinogram, [], 1), 1);
filteredProj = ft_R .* myFilter;
filteredProj = ifftshift(filteredProj, 1);
sinogram = real(ifft(filteredProj,[],1));
sinogram = sinogram';
end