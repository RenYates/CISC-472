% Simulate CT scanning and reconstruction

nTheta = 200; % number of discrete angles in sinogram

% Spatial filters to apply to backprojections

RamLakFilter = [0, -1/(9*pi*pi), 0, -1/(pi*pi), 0.25, -1/(pi*pi), 0, -1/(9*pi*pi), 0];
RamLakFilter = RamLakFilter ./ sum(RamLakFilter);

SheppLoganFilter = arrayfun( @(i) (-2/(pi*pi*(4*i*i-1))), [-4, -3, -2, -1, 0, 1, 2, 3, 4] );
SheppLoganFilter = SheppLoganFilter ./ sum(SheppLoganFilter);

% Tiled output window

t = tiledlayout(4,2); % 4x2 images
t.TileSpacing = 'compact';
t.Padding = 'compact';
 
% Read image.  Map each greyscale value to [0,1]

% img = im2double( rgb2gray( imread( 'skull.png' ) ) );
img = im2double( rgb2gray( imread( 'SheppLogan_Phantom.png' ) ) );

nexttile; imshow(1); nexttile; imshow(img); title( 'Original image' );

% Build sinogram

S = parallelSinogram( img, nTheta );

sinoFactor = 1 / max(max(S));  % Scale sinograms so that max value = 1 when displaying.
                               % Use same factor for all sinograms.

nexttile; imshow(S .* sinoFactor ); title( 'Unfiltered sinogram' );

% Backproject unfiltered sinogram

BP1 = backproject( S, size(img,1) );

nexttile; imshow(BP1); title( 'Unfiltered backprojection' );

% Apply Ram-Lak filter to sinogram

S2 = zeros( size(S) );

for r = 1:size(S,1)
    S2(r,:) = applyFilter( RamLakFilter, S(r,:) );
end

nexttile; imshow(S2 .* sinoFactor); title( 'Ram-Lak-filtered sinogram' );

% Backproject Ram-Lak-filtered sinogram

BP2 = backproject( S2, size(img,1) );

nexttile; imshow(BP2); title( 'Ram-Lak-filtered backprojection' );

% Apply Shepp-Logan filter to sinogram

S3 = zeros( size(S) );

for r = 1:size(S,1)
    S3(r,:) = applyFilter( SheppLoganFilter, S(r,:) );
end

nexttile; imshow(S3 .* sinoFactor); title( 'Shepp-Logan-filtered sinogram' );

% Backproject Shepp-Logan-filtered sinogram

BP3 = backproject( S3, size(img,1) );

nexttile; imshow(BP3); title( 'Shepp-Logan-filtered backprojection' );

% Apply iterative reconstructions

numIterations = 50; % = number of individual *rays* to correct
numIntermediateImages = 0; % can show images during reconstruction, but need space in tiledlayout() at top

BP4 = iterativeReconstruction( S, size(img,1), numIterations, numIntermediateImages );

nexttile(1); imshow(BP4); title( sprintf( '%d iterations', numIterations ) );


% Return the parallel-ray sinogram of image 'img' with 'nTheta' angles
% and a number of samples at each angle each to the number of image rows.
%
% Each entry in the sinogram is the length-weighted sum of attenuation
% values (from [0,1]) along a ray.  So the sum can be greater than 1.
%
% The image should be square and the subject should not go outside the
% largest circle inside the square, as is the case with CT scans.

function S = parallelSinogram( img, nTheta )

    % check that image is square
    
    dims = size(img);
    if dims(1) ~= dims(2)
        error( 'Image dimensions ' + dims + ' are not equal.  Image must be square.' );
    end

    % set up sinogram
    
    nRho = dims(1);

    S = zeros( nTheta, nRho );
    
    % calculate each row of sinogram
    
    % Determine the amount of rotation between each theta
    rotationIncrement = 180/nTheta;
    
    % For each theta value
    for i = 1:(nTheta)
        
       theta = rotationIncrement*(i-1);
       
       %Rotate the image by theta
       imRotated = imrotate(img,theta, 'bilinear','crop');
       
       %Sum the values of the image along the first dimension
       rho = sum(imRotated);
       
       %Add the new rho to the ith row of the sinogram
       S(i,:) = rho;
                     
    end

end


% Apply the Ram-Lak filter to a sinogram row.
%
% filteredRow(i) = sum_j filter(j) * row(i-j+offset)
%
% 'row' is zero outside its range.

function filteredRow = applyFilter( filter, row )
    
    offset = (length(filter)+1)/2; % offset to centre of filter
    
    filteredRow = zeros( 1, length(row) );

    %For each column in row
    for i=1:length(row)
        
        tempSum = 0;
        
        %For each value in filter
        for j=1:length(filter)
            
            %Verify that row index is within row dimensions
            if (i-j+offset) >= 1 && (i-j+offset) <= length(row)
  
                %Sum applied filter values on row
                tempSum = tempSum + (filter(j)*row(i-j+offset));
                
            end
        end
        
        %Update ith position in filteredRow with filtered value
        filteredRow(i) = tempSum;
        
    end
    
end


% Return a 'dim' x 'dim' image that is the backprojection of sinogram 'S'.
%
% Each sinogram value is normalized by dividing by the number of pixels across 
% which it is backprojected, and divided by the number of angles that backprojections
% are made from.  The number of pixels here should be the same as the number of 
% pixels used when computing the sinogram.

function BP = backproject( S, dim )
    
    % set up backprojection image

    BP = zeros( dim, dim );

    % for each angle, project sinogram rays across backprojection

    nTheta = size(S,1);

    % Determine the amount of rotation between each theta
    rotationIncrement = 180/nTheta;
    
    %Go through each theta value
    for i = 1:(nTheta)
                
       rho = S(i,:);
       
       %Smear rho values on temporary image
       tempIm = repmat(rho, dim,1);
       
       %Normalize rho values against the number of pixels
       normTempIm = tempIm/dim;
       
       %Calculate rotation for current theta
       theta = rotationIncrement*(i-1);   
       
       %Rotate the image by negative theta
       imRotated = imrotate(normTempIm,-1*theta, 'bilinear','crop');
       
       %Add smeared values to backprojection collector array
       BP = BP + imRotated;
                     
    end
    
    %Divide BP values by number of rotations
    BP = BP/nTheta;
    
    %Normalize BP image for better visualization
    BP = uint8(255*mat2gray(BP));
    
end


% Perform iterative reconstruction

function IR = iterativeReconstruction( S, dim, numIterations, numImagesOutput )
    
    % set up reconstructed image
    
    IR = zeros( dim, dim );

    nTheta = size(S,1);
    
    iterationsComplete = 0;
    iterationsPerImage = idivide( numIterations, int32(numImagesOutput) );
        
    while 1
        for i = randperm(nTheta) % do angles in random order
            
            if iterationsComplete >= numIterations
                return;
            end
            
            %Determine current theta value
            theta = (180/nTheta)*i;
            
            %Rotate IR by theta
            IRRotated = imrotate(IR, theta, 'bilinear','crop');
            
            %Calculate current rho estimate at theta
            rhoGuess = sum(IRRotated);
            
            %Determine actual value of rho at theta
            rhoActual = S(i,:);
            
            %Calculate difference
            rhoDiff = rhoGuess - rhoActual;
            
            %Divide difference by number of pixels in dim
            rhoDiffNorm = rhoDiff/dim;
            
            %Apply rho difference to rotated iterative reconstruction image
            for j = 1:dim
                IRRotated(:,j) = IRRotated(:,j) - rhoDiffNorm(j);
            end
            
            %Rotate image back and set equal to IR
            IR = imrotate(IRRotated, -1*theta, 'bilinear','crop');
        
            iterationsComplete = iterationsComplete + 1;
            
            if mod( iterationsComplete, iterationsPerImage ) == 0
                nexttile; imshow(IR); title( sprintf( '%d iterations', iterationsComplete ) );
            end
        end
    end
end


