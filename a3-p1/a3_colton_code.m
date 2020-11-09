% Perform ICP to align collected points with a model

modelFilename = 'femur.stl';            % model
collectedFilename = 'knee1.csv';        % YOUR OWN collected points
stylusID = '8700340';                   % ID in collected points file

% ICP parameters

maxRMSError = 5.0;             % ICP error (in mm)
maxIterations = 200;           % ICP iterations
noImprovementFraction = 0.001; % quit iterating when RMSE improves by less than this fraction
numAttempts = 5;               % number of attempts with ICP
randomICP = false;

% PROVIDE YOUR OWN CALIBRATION HERE FROM ASSIGNMENT 2

% calculated using pivot_calibration_1 file from assignment 2
stylusTip = [-17.4394 -0.998029 -159.076];  % stylus tip in stylus coordinate system

% Get the stylus-collected points

pts = read_NDI_data( collectedFilename, stylusID, stylusTip );
disp( sprintf( 'read %d stylus points', length(pts)) );

observedPts = find_observed_points( pts );
n = size(observedPts,1);

disp( sprintf( 'found %d observed stylus points', length(observedPts)) );

% read the model and get its unique vertices

[modelFaces, modelPts] = stlread( modelFilename );
modelPts = unique( modelPts, 'rows' ); % remove duplicates

disp( sprintf( 'read %d model points', length(modelPts)) );

% Build KD-tree from vertices

kdTree = KDTreeSearcher( modelPts, 'BucketSize', 10 );

%Specify the range within which the center of mass can exist
xRange = [-150,50];
yRange = [-150,-100];
zRange = [-350,-250];

%Calculate the mean of the observed points
observedMean = mean(observedPts);

%Collect the best accumRot and accumTrans seen so far
bestAccumRot = eye(3,3);
bestAccumTrans = [0 0 0];
bestRMSE = 9999;

for i = 1:numAttempts
    
    % Find ICP transformation
    
    rmsError = 9999;    % current error
    iter = 1;           % current iteration

    % Accumulate a transformation over multiple iterations
    %
    % This transformation will transform OBSERVED points onto MODEL points.

    accumRot   = eye(3,3);  % accumulated rotation 
    accumTrans = [0 0 0];   % accumulated translation
    prevRMSE = 0;
    
    if randomICP == true
        rotAngleX = abs((0-360)*rand(1));
        rotAngleY = abs((0-360)*rand(1));
        rotAngleZ = abs((0-360)*rand(1));
    
        %Set a random initial rotation
        rotMatrixX = [1 0 0; 0 cosd(rotAngleX) -1*sind(rotAngleX); 0 sind(rotAngleX) cosd(rotAngleX)];
        rotMatrixY = [cosd(rotAngleY) 0 sind(rotAngleY); 0 1 0; -1*sind(rotAngleY) 0 cosd(rotAngleY)];
        rotMatrixZ = [cosd(rotAngleZ) -1*sind(rotAngleZ) 0; sind(rotAngleZ) cosd(rotAngleZ) 0; 0 0 1];
        initRot = rotMatrixZ*rotMatrixY*rotMatrixX;
    
        %Set a random initial translation
        randX = (xRange(2)-xRange(1))*rand(1)+xRange(1);
        randY = (yRange(2)-yRange(1))*rand(1)+yRange(1);
        randZ = (zRange(2)-zRange(1))*rand(1)+zRange(1);
        randStartLocation = [randX, randY, randZ];
    
        %Determine the vector from observedMean to randStartLocation
        initTrans = randStartLocation - (initRot*(observedMean.')).';
    
        %Apply this translation to xPts
        initialPts = (initRot*(observedPts.') + initTrans.').';
        xPts = initialPts; % transformed points (initially = observed points)
    
        accumRot = initRot * accumRot;
        accumTrans = ((initRot * accumTrans.') + initTrans.').';
    else
        xPts = observedPts;
    end

    while rmsError > maxRMSError & iter < maxIterations & abs(rmsError - prevRMSE) > noImprovementFraction*prevRMSE

      prevRMSE = rmsError;

      % Find the nearest model points using kD-tree

      indices = knnsearch( kdTree, xPts );
      closestPts = modelPts(indices,:);

      draw_all( modelPts, xPts, closestPts );

      % calculate centre of point clusters
      xPtsMean = mean(xPts);
      closestPtsMean = mean(closestPts);

      xPtsZeroMean = xPts - xPtsMean;
      closestPtsZeroMean = closestPts - closestPtsMean;

      % implement procrustes

      % calculate covariance matrix by multiplying col vectors xPts with row
      % vectors in closestPts (use zero mean for svd)
      CovarianceMatrix = xPtsZeroMean.' * closestPtsZeroMean;

      % calculate singular value decomposition
      % U is the orthogonal basis for columns in CovarianceMatrix
      % V is the orthogonal basis for rows in CovarianceMatrix
      [U,~,V] = svd(CovarianceMatrix);

      R = V * U.'; % rotation matrix from xPts to closestPts

      t = closestPtsMean.' - (R * xPtsMean.'); % translation vector from xPts to closestPts

      % accumulate rotation and translations
      accumRot = R * accumRot;
      accumTrans = ((R * accumTrans.') + t).';

      % update xPts by rotating and translating them
      xPts = (accumRot*(observedPts.') + accumTrans.').';

      % calculate rmsError
      rmsError = sqrt(mean(norm(closestPts - xPts).^2));

      % Update the display
      disp( sprintf( '%2d: RMSE = %.2f', iter, rmsError ) );

      iter = iter + 1;
    end
    
    if rmsError < bestRMSE
        bestAccumRot = accumRot;
        bestAccumTrans = accumTrans;
        bestRMSE = rmsError;
        bestClosestPts = closestPts;
    end

% Final rendering and report

draw_all( modelPts, xPts, closestPts );

if abs(rmsError - prevRMSE) <= noImprovementFraction*prevRMSE
  disp( sprintf( 'Stopped due to insufficient (%.2f%%) improvement.', abs(rmsError-prevRMSE)/prevRMSE*100 ) );
elseif iter >= maxIterations
  disp( sprintf( 'Stopped after limit of %d iterations was reached.', maxIterations ) );
else	
  disp( 'Stopped with low RMS error.' );
end

% Show resulting rotation and translation

%accumRot
%accumTrans

end

% draw best result to save as picture example
bestAccumRot
bestAccumTrans
bestRMSE

bestPts = (bestAccumRot*(observedPts.') + bestAccumTrans.').';

draw_all( modelPts, bestPts, bestClosestPts );
disp("best result has been drawn");

% ---- Done ----



% In a stream of 3D points, find those that are in approximately the same
% position for a "long time".

function ptsOut = find_observed_points( ptsIn )

  minRestingPts = 20;   % need at least this many observed resting points
  maxRestingDist = 2.0; % max distance (mm) between observed resting points
  
  ptsOut = [];
  pt0 = [0 0 0];
  count = 0;
  
  for k = 1:length(ptsIn) % Note: this will discard the last set of resting points
                          % when it's not followed by a non-resting point
    pt = ptsIn(k,:);
    dist = norm( pt-pt0 );
    if dist <= maxRestingDist % add a point
      count = count+1;  
    elseif count < minRestingPts % end of resting points ... too few
      pt0 = pt;
      count = 0;
    else % end of resting points ... enough
      ptsOut = [ ptsOut; median( ptsIn( k-count:k-1, : ) ) ];
      count = 0;
    end
  end
end


% Draw the model points, the current (transformed) observed points,
% and lines from the observed points to their closest model points.

function draw_all( modelPts, observedPts, closestPts )
    
  f = figure(1);
  clf(f);
  
  axis('vis3d');
    
  % Draw the model points
  
  scatter3( modelPts(:,1), modelPts(:,2), modelPts(:,3), 1, [0.8, 0.8, 0.8] ); % grey

  % Draw the transformed observed points
  
  hold on;
  
  scatter3( observedPts(:,1), observedPts(:,2), observedPts(:,3), 30, [0.8, 0.3, 0.1] ); % reddish
  
  % Draw lines between correspondings points
  
  for i = 1:length(observedPts)
    plot3( [observedPts(i,1), closestPts(i,1)], [observedPts(i,2), closestPts(i,2)], [observedPts(i,3), closestPts(i,3)] ); 
  end
  
  drawnow;
  hold off;
end


% Read a CSV file from the NDI tracker software and extract the
% position and orientation columns of the named marker.
%
% pos    = n x 3 of translations
% orient = n x 4 of quaternion orientations

function pts = read_NDI_data( dataFile, markerID, localVector )

  t = readtable( dataFile, 'PreserveVariableNames', true );

  % Find the column 

  colIndex = find(contains( t.Properties.VariableNames, markerID ));

  if colIndex == []
    disp( sprintf( "In %s: Could not find a column header containing '%s'.", dataFile, markerID ) );
    exit;
  end

  % From the rows with state == 'OK' (state is at +3 offset from the ID
  % column), extract the Tx, Ty, Tz columns.

  status = t{:,colIndex+3};
  n = size( t( strcmp(status,'OK'), 1 ), 1 );

  pts = zeros( n, 3 );
  k = 1;

  for i = 1:size(t,1)
    if strcmp( t{i,colIndex+3}, 'OK' )

      % Coerce the columns to 'double', since 'readtable' sometimes
      % records numbers as strings.  MATLAB BUG!

      % Extract pose's rotation as a quaternion.  This is in columns
      % offset by +4, +5, +6, +7 from the ID column.

      quat = zeros( 1, 4 );
      for j=1:4
        if iscell( t{i,colIndex+3+j} )
          quat(j) = str2double( t{i,colIndex+3+j}{1} );
        else
          quat(j) = t{i,colIndex+3+j};
        end
      end

      % Extract pose's translation as a vector.  This is in columns
      % offset by +8, +9, +10 from the ID column.

      pos = zeros( 1, 3 );
      for j=1:3
        if iscell( t{i,colIndex+7+j} )
          pos(j) = str2double( t{i,colIndex+7+j}{1} );
        else
          pos(j) = t{i,colIndex+7+j};
        end
      end
      
      pts(k,:) = (quaternion_to_matrix( quat ) * localVector')' + pos;

      k = k + 1;
    end
  end
end


% From https://www.mathworks.com/matlabcentral/fileexchange/35475-quaternions
%
% Convert quaternion to 3x3 matrix.

function R = quaternion_to_matrix( Qrotation )

  w = Qrotation( 1 );
  x = Qrotation( 2 );
  y = Qrotation( 3 );
  z = Qrotation( 4 );

  Rxx = 1 - 2*(y^2 + z^2);
  Rxy = 2*(x*y - z*w);
  Rxz = 2*(x*z + y*w);
  Ryx = 2*(x*y + z*w);
  Ryy = 1 - 2*(x^2 + z^2);
  Ryz = 2*(y*z - x*w );
  Rzx = 2*(x*z - y*w );
  Rzy = 2*(y*z + x*w );
  Rzz = 1 - 2 *(x^2 + y^2);

  R = [ Rxx,    Rxy,    Rxz;
        Ryx,    Ryy,    Ryz;
        Rzx,    Rzy,    Rzz  ];
end
