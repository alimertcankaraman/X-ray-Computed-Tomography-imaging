clear; clc;

%% INPUTS

imageName = input('Enter the name of the image: ','s');
img = load(imageName);

[M, N] = size(img.(imageName));

numberOfAngles = str2double(input('Enter the number of projections: ', 's'));
numberOfBeams = str2double(input('Enter the number of sampling points in each projection: ', 's'));

stepSize = (180 - 0) / numberOfAngles;
thetaValues = linspace(0, (180 - stepSize), numberOfAngles);
tValues = linspace((-M*sqrt(2)/2), (M*sqrt(2)/2), numberOfBeams);

xCoor = floor(-M/2):ceil(M/2);
yCoor = floor(-N/2):ceil(N/2);

pt = zeros(numberOfAngles,numberOfBeams);
distances_teta = zeros(numberOfBeams, (length(xCoor)*2-1), numberOfAngles);
rows_teta = zeros(numberOfBeams, (length(xCoor)*2-1), numberOfAngles);
cols_teta = zeros(numberOfBeams, (length(xCoor)*2-1), numberOfAngles);

%% FORWARD PROJECTION

for teta = 1:numberOfAngles 

    distances = zeros(numberOfBeams, (length(xCoor)*2-1));
    rows = zeros(numberOfBeams, (length(xCoor)*2-1));
    cols = zeros(numberOfBeams, (length(yCoor)*2-1));
    
    if(0 == thetaValues(teta))
        %Find the midpoints and length of the line segments ignoring the irrelevant points
        for t =  1:numberOfBeams
            for i = 1:(length(yCoor)-1)
                %Calculate length and address of the line segment if it intersects the image
                if(M/2 >= abs(tValues(t)))  
                    distances(t,i) = yCoor(i+1) - yCoor(i); 
                    rows(t,i) = M/2 - yCoor(i);      
                    cols(t,i) = M/2 + ceil(tValues(t));
                end
                %Calculation of projection function
                if((0 < rows(t,i)) && (0 < cols(t,i)))
                    pt(teta,t) = pt(teta,t) + img.(imageName)(rows(t,i), cols(t,i)) * distances(t,i);
                end
            end
        end
    elseif(90 == thetaValues(teta))
        %Find the midpoints and length of the line segments ignoring the irrelevant points
        for t =  1:numberOfBeams
            for i = 1:(length(xCoor)-1)
                %Calculate length and address of the line segment if it intersects the image
                if(M/2 >= abs(tValues(t)))  
                    distances(t,i) = xCoor(i+1) - xCoor(i); 
                    cols(t,i) = M/2 - xCoor(i);                     
                    rows(t,i) = M/2 - floor(tValues(t));
                end
                %Calculation of projection function
                if((0 < rows(t,i)) && (0 < cols(t,i)))
                    pt(teta,t) = pt(teta,t) + img.(imageName)(rows(t,i), cols(t,i)) * distances(t,i);
                end
            end
        end
    else
        xValues = zeros(1, length(xCoor));
        yValues = zeros(1, length(yCoor));
        midXPoints = zeros(1, length(xCoor)-1);
        midYPoints = zeros(1, length(yCoor)-1);
        tempx = zeros(1,length(xCoor)*2);
        tempy = zeros(1,length(xCoor)*2);
        xPoints = zeros(1,length(xCoor)*2);
        yPoints = zeros(1,length(xCoor)*2);
        combined = zeros(length(xCoor)*2,2);
        points = zeros(length(xCoor)*2,2);

        for t =  1:numberOfBeams
            for i = 1:(length(xCoor))
                %Find the x and y intersection points
                xValues(i) = (tValues(t) - yCoor(i) * sind(thetaValues(teta))) / cosd(thetaValues(teta));
                yValues(i) = (tValues(t) - xCoor(i) * cosd(thetaValues(teta))) / sind(thetaValues(teta));
            end

            tempx = [xCoor xValues];
            tempy = [yValues yCoor];
            xPoints = tempx;
            yPoints = tempy;

            %Check if the x and y points are on the image otherwise, remove from the matrix
            xPoints(abs(tempx) > M/2) = [];
            xPoints(abs(tempy) > M/2) = [];
            yPoints(abs(tempx) > M/2) = [];
            yPoints(abs(tempy) > M/2) = [];

            %1st column is x points and 2nd column is y points
            combined = [xPoints' yPoints'];

            %Sort the relevant points
            points = unique(round(sortrows(combined, 1), 5), 'rows');
        
            if(1 < (size(points,1)))
                for i = 1:(size(points,1)-1)
                    %Calculation of distance, row and column data
                    distances(t,i) = sqrt( (abs(points(i+1,1)) - abs(points(i,1)))^2 + ...
                                           (abs(points(i+1,2)) - abs(points(i,2)))^2 );
                    midXPoints(i) = (points(i,1) + points(i+1,1)) / 2;
                    midYPoints(i) = (points(i,2) + points(i+1,2)) / 2;
                    rows(t,i) = M/2 - floor(midYPoints(i));
                    cols(t,i) = M/2 + ceil(midXPoints(i)); 
                    
                    %Calculation of projection function
                    if((0 < rows(t,i)) && (0 < cols(t,i)))
                        pt(teta,t) = pt(teta,t) + img.(imageName)(rows(t,i), cols(t,i)) * distances(t,i);
                    end 
                end
            end
        end
    end

    %Store each distance, row and column data for each projection
    distances_teta(:,:,teta) = distances; 
    rows_teta(:,:,teta) = rows; 
    cols_teta(:,:,teta) = cols; 
end

%% VISUALIZATION
% Comparison with the radon() function.

[y,xp] = radon(img.(imageName),thetaValues);

for teta = 1:length(thetaValues)
    ptNormalized(teta,:) = pt(teta,:) / max(pt(teta,:));
    yNormalized(:,teta) = y(:, teta) / max(y(:,teta));
end

projectionNumber = round(numberOfAngles/2);
figure;
g(1) = plot(pt(projectionNumber,:), "-r");
hold on
g(2) = plot(y(:,projectionNumber),"-b");
grid on
xlabel("t Values");
ylabel("p(t)");
legend(g, 'myOutput', 'radon()')
title("theta:", thetaValues(projectionNumber));
%}

%% FILTER

projectionsFiltered = fRamlak(pt);

%% BACKPROJECTION

image = zeros(M,N);
imageFiltered = zeros(M,N);
distances = zeros(numberOfBeams, (length(xCoor)-1));
rows = zeros(numberOfBeams, (length(xCoor)-1));
cols = zeros(numberOfBeams, (length(yCoor)-1));
ptt = zeros(numberOfBeams);
ptFiltered = zeros(numberOfBeams);

for teta = 1:numberOfAngles
    distances = distances_teta(:,:,teta);
    rows = rows_teta(:,:,teta);
    cols = cols_teta(:,:,teta);
    ptt = pt(teta, :);
    ptFiltered = projectionsFiltered(teta, :);

    for t = 1:numberOfBeams
        for i = 1:(length(xCoor)-1)
            if((0 < rows(t,i)) && (0 < cols(t,i)))
                image(rows(t,i), cols(t,i)) = image(rows(t,i), cols(t,i)) + ptt(t) * distances(t,i);
                imageFiltered(rows(t,i), cols(t,i)) = imageFiltered(rows(t,i), cols(t,i)) + ptFiltered(t) * distances(t,i);
            end
        end
    end
end

%% VISUALIZATION

normalized = img.(imageName) / max(img.(imageName)(:));
imageNormalized = image / max(image(:));
imageFilteredNormalized = imageFiltered / max(imageFiltered(:));

%[R,xp] = radon(normalized,thetaValues);
I = iradon(pt', thetaValues','linear','none');
I = I / max(I(:));
I2 = iradon(pt', thetaValues');
I2 = I2 / max(I2(:));

figure;
imshow(normalized);
title('Original Image')
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
