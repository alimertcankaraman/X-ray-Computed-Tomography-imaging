function sinogram = filterRamlak(sinogram)
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