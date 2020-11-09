% Suppress warnings
warning('off','all');

% Use sphere fitting to calibrate pointer 8700340

stylusID = '8700340';
markerID = '8700339';
% uncomment line to read data file individually
%dataFile = 'pivot_calibration_0.csv';
dataFile = 'pivot_calibration_1.csv';
%dataFile = 'pivot_calibration_2.csv';

setupDrawing();

% Read the raw data into 'pos' (a translation 3-vector) and 'orient'
% (a quaternion 4-vector).

[pos, orient] = read_NDI_data( dataFile, stylusID );

% fit the sphere to get centre c and radius r

%[c, r] = fitSphere( pos );

% uncomment lines 25-29 for RANSAC sphere fit
[c, r, inlierInd] = fitSphereWithRANSAC(pos, orient);

% update pos and orient to only include inlier points
%pos = pos(inlierInd,:);
%orient = orient(inlierInd,:);

% Show the fit
drawCoordSystems( pos, orient );
drawSphere( c, r );

% Transform c into the coordinate system of each pose
numPoses = size(pos,1);
allC = zeros(numPoses, 3); % row vectors

for i=1:numPoses
    rotation = quaternion_to_matrix(orient(i,:)); % stylus to world
    t = pos(i,:);
    allC(i,:) = rotation.' * (c.' - t).'; % stylus tip in pose coord systems
end

% Find the average transformed c, which should be the same in all of
% the stylus coordinate systems.  Also find the standard deviation.

averageC = mean(allC,1);
stdevC = std(allC);

% Report the results
%
% 'c_average' is the average tip position in the stylus coordinate
% system.  'c_stdev' is its standard deviation in the stylus
% coordinate system.

c_average = averageC;
c_stdev = stdevC;


disp( sprintf( 'tip position in stylus CS: (%g, %g, %g)', c_average(1), c_average(2), c_average(3) ) );
disp( sprintf( 'tip stdev in stylus CS:    (%g, %g, %g)', c_stdev(1),   c_stdev(2),   c_stdev(3) ) );


% Show vectors to tips in world coordinate system
%
% This is for debugging, so that you can see that the vector touches
% the same pivot point from all stylus coordinate systems.

drawLocalVectorInCoordSystems( c_average, pos, orient );

% Show tip points in global system, along with 95% confidence
% interval as an ellipsoid.
%
% 'c_world' are the tip points in the world coordinate system.
% They should all be very near the pivot point.

numPoses = size(pos,1);
c_world = zeros(numPoses, 3); % row vectors

for i=1:numPoses
    rotation = quaternion_to_matrix(orient(i,:)); % stylus to world
    t = pos(i,:);
    c_world(i,:) = rotation*(c_average.') + t.'; % stylus tip in pose coord systems
end

c_world_stdev = std(c_world);
c_world_average = mean(c_world);

disp( sprintf( 'tip position in world CS: (%g, %g, %g)', c_world_average(1), c_world_average(2), c_world_average(3) ) );
disp( sprintf( 'tip stdev in world CS:    (%g, %g, %g)', c_world_stdev(1),   c_world_stdev(2),   c_world_stdev(3) ) );

drawPointsWithEllipsoid( c_world, c_world_stdev ); % pass stdev of c_world

% ---------------- END OF MAIN CODE ----------------

% Fit a sphere to a set of positions
%
% See http://watkins.cs.queensu.ca/~jstewart/472/notes/08-matrix-applications/08-matrix-applications.html

function [c, r] = fitSphere( pos )
    row = size(pos,1);
    b = zeros(row,1);
    A = zeros(row,4);
    % Calculate A and b matrices in Ax=b equation using point positions
    for index=1:row
        b(index,:) = dot(pos(index,:),pos(index,:));
        A(index,:) = [2*pos(index,:), 1];
    end
    % Solve for x in Ax = b equation (to find center and radius of sphere)
    x = pinv(A)*b;
    c = x(1:3); % center point is first three values in 4x1 vector x
    r = sqrt(x(4)+dot(c,c)); % solve for r using center point
end
  

% Fit a sphere to a set of positions using RANSAC.
%
% ALSO RETURN THE INDICES OF THE BEST INLIERS.  THE CALLING CODE
% SHOULD RESTRICT ITSELF TO THOSE INLIERS.
%
% See https://en.wikipedia.org/wiki/Random_sample_consensus

function [c, r, bestInlierIndices] = fitSphereWithRANSAC( pos, orient )

    % acceptance values for RANSAC estimation to exit early or update
    ppAcceptance = 0.90;
    pointPercent = 0;
    radiusPercent = 0.10;
    bestInlierIndices = [];
    
    rows = size(pos,1);
    
    for i=0:200 % set number of iterations to prevent infinite loop
        inlierIndices = [];
        % randomly generate numbers
        randNums = randperm(rows);
        % grab the first 4 column values to randomly select points from pos
        points = pos(randNums(1:4),:);
        [c, r] = fitSphere(points);
        %inlier points must be within a percentage of radius
        radMax = r + r*radiusPercent;
        radMin = r - r*radiusPercent;
        
        % iterate through all points in pos
        for index=1:rows
            % calculate the distance between the center of the sphere to
            % the points
            dpoint = norm(c - pos(index,:).');
            % if the point is within % range around radius, set as inlier
            if dpoint < radMax && dpoint > radMin
                inlierIndices = [inlierIndices; index];
            end
        end
        
        % calculate percentage of points included as inliers
        numInliers = size(inlierIndices,1);
        tempPercent = numInliers / rows;
        
        % if the percentage of points is > the threshold, update values
        % and exit loop
        if tempPercent > ppAcceptance
            disp("loop ended early")
            pointPercent = tempPercent;
            bestInlierIndices = inlierIndices;
            drawCoordSystems( pos(inlierIndices,:), orient(inlierIndices,:) );
           break;
        % if the percentage of points is > than current best, update values
        elseif tempPercent > pointPercent
            disp("inliers updated")
            pointPercent = tempPercent;
            bestInlierIndices = inlierIndices;
        end
    end
    inlierPoints = pos(bestInlierIndices,:);
    
    % refit a sphere to all the inlier points
    [c, r] = fitSphere(inlierPoints);
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


% Read a CSV file from the NDI tracker software and extract the
% position and orientation columns of the named marker.
%
% pos    = n x 3 of translations
% orient = n x 4 of quaternion orientations

function [pos, quat] = read_NDI_data( dataFile, markerID )

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

  pos = zeros( n, 3 );
  quat = zeros( n, 4 );
  k = 1;

  for i = 1:size(t,1)
    if strcmp( t{i,colIndex+3}, 'OK' )

      % Coerce the columns to 'double', since 'readtable' sometimes
      % records numbers as strings.  MATLAB BUG!

      % Extract pose's rotation as a quaternion.  This is in columns
      % offset by +4, +5, +6, +7 from the ID column.

      for j=1:4
	if iscell( t{i,colIndex+3+j} )
	  quat(k,j) = str2double( t{i,colIndex+3+j}{1} );
	else
	  quat(k,j) = t{i,colIndex+3+j};
	end
      end

      % Extract pose's translation as a vector.  This is in columns
      % offset by +8, +9, +10 from the ID column.

      for j=1:3
	if iscell( t{i,colIndex+7+j} )
	  pos(k,j) = str2double( t{i,colIndex+7+j}{1} );
	else
	  pos(k,j) = t{i,colIndex+7+j};
	end
      end

      k = k + 1;
    end
  end
  
  disp( sprintf( '%d points collected', size(pos,1) ) );
end



% Set up the drawing

function setupDrawing() 

    f = figure(1);
    clf(f);
    view(3);
    daspect( [1 1 1] );
    pbaspect manual;
    hold on;
end


% Draw a set of coordinate systems
%
% pos and orient store their vectors in rows.

function drawCoordSystems( pos, orient )
    
    colours = [ 'r' 'g' 'b' ];

    scale = 0.005 * norm(max(pos) - min(pos));

    for i=1:size(pos,1)
        m = quaternion_to_matrix( orient(i,:) );
        t = pos(i,:);
        for j=1:3
            head = t + m(:,j)' .* scale;
            plot3( [t(1) head(1)], [t(2) head(2)], [t(3) head(3)], colours(j) );
        end
    end
end


% Draw a sphere

function drawSphere( c, r )
    
    [x,y,z] = sphere;
    surf( x*r+c(1), y*r+c(2), z*r+c(3), 'FaceAlpha', 0.05, 'FaceColor', [0.6 0.3 0.3] );
end


% Draw a local vectors in different coordinate systems
%
% v is a row vector.  pos and orient store their vectors in rows, too.

function drawLocalVectorInCoordSystems( v, pos, orient )
    
    for i=1:size(pos,1)
        m = quaternion_to_matrix( orient(i,:) );
        t = pos(i,:);
        for j=1:3
            head = t + (m * v')';
            plot3( [t(1) head(1)], [t(2) head(2)], [t(3) head(3)] );
        end
    end
end


% Draw points with a 95% CI ellipsoid.
%
% Use matlab's 'ellipsoid' and 'surf' functions.
%

function drawPointsWithEllipsoid( points, stdev )
    figure;
    hold on;
    
    % Professor Stewart's code for aligning ellipsoid on point cloud axis
    [X,Y,Z] = ellipsoid(0,0,0,1,1,1);

    [eigvec,eigval] = eig(cov(points));

    XYZ = [X(:),Y(:),Z(:)] * (1.96*sqrt(eigval)) * eigvec';

    mu = mean( points );

    X(:) = XYZ(:,1)+mu(1);
    Y(:) = XYZ(:,2)+mu(2);
    Z(:) = XYZ(:,3)+mu(3);

    surf( X, Y, Z, 'FaceAlpha', 0.1);

    mn = min(points);
    mx = max(points);
    axis( [mn(1),mx(1),mn(2),mx(2),mn(3),mx(3)] );

    %This draws the vector along the stylus axis:

    tip = mean( points); 
    arrow3( [0, 0, 0], [tip(1), tip(2), tip(3)], 'f', 1 );
    
    % plot best fit vectors as points
    scatter3(points(:,1), points(:,2), points(:,3),5, "filled");
end
