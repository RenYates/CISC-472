Group Members:

Colton Barr, 10161906, 14cb49
Lauren Yates, 10195969, 15ly1

------------------------------------------------------

#7 - Answers

The main differences between our own pivot data and the pivot data provided were the pivot location relative to the tracker, the presence of tracking errors, and the number of outliers. Since both of our pivot calibrations were captured with the same pivot location relative to the tracker, the calculated spheres and center points were almost identical between our own trials. The provided data was performed with a different pivot point, however, and produced a different center point in tracker coordinates. Within the data we captured, we also recorded frames where the stylus reflectors were obscured from the tracker’s view. These appeared as rows in the .csv with the values “MISSING” in all dimensions, and were removed prior to use to prevent errors related to reading in invalid data. We also captured one trial with a significantly smaller number of outliers than the other two files; this occurred because tracking was started and stopped with the probe in the correct pivot location.

The core purpose of RANSAC is to deal with outliers in a dataset, and we found that in the 2 datasets with large outliers RANSAC consistently outperformed regular sphere fitting. For both the dataset provided and our second pivot calibration, the standard deviation of the tip position in world coordinates was lower using RANSAC than using simple sphere fitting. This occurred because the inlier points detected by RANSAC excluded the major outliers, which generally corresponded to frames where the pointer was not on the same pivot point. Fitting the sphere to these inliers was virtually equivalent to fitting the sphere exclusively to frames in which the probe was within the pivot dimple.

In the case of the first pivot calibration with minimal outliers, we found that RANSAC performed identically to regular sphere fitting. This occurred because all points were detected as inliers by RANSAC and used for the subsequent sphere fitting step, resulting in an identical fit to the non-RANSAC sphere fitting approach. Since RANSAC is principally designed for outlier removal, it is intuitive that using RANSAC for data free of outliers would not afford any discernable advantage. Clearly in the case of very clean data RANSAC adds an unnecessary step and may serve to slow down the process of sphere fitting.

The main recommendation we would have related to pivot calibration is to limit the number of outliers generated during data capture. Outliers can be limited by ensuring that pivoting is performed in a small, well defined divot on a surface that is rigid relative to the tracker. That is to say, the location of the pivot point should not change relative to the location of the tracker over the course of calibration. The capturing of points should also start and end with the stylus already in the calibration divot, and care should be taken to ensure the probe does not leave the divot during calibration. If a large volume of outliers are captured during calibration, the user should recapture the points to avoid skewing the data and increasing the standard deviation. Since RANSAC does not negatively influence a calibration with a low number of outliers, it should be used regardless of the quality of the points captured. This will ensure that if outliers are present they are excluded as is deemed appropriate.

Since the 4 retroreflectors on the stylus are coplanar, we would expect the error in tracking the stylus to be highest perpendicular to the plane of the retroreflectors. This is because a change in position perpendicular to the plane of the retroreflectors will appear only as a subtle change in the relative positions on the camera sensors. A change in position in the plane of the retroreflectors, however, manifests as a much larger change in the position of the retroreflectors in the camera sensors. Since reliable detection of the pointer relies on maintaining line of sight between the reflectors and the camera, the tool must be manipulated with the plane of the reflectors relatively coplanar with the plane of the camera sensors. This means that we should expect to see the greatest error in localization in the axis perpendicular to the plane of the reflectors, and this axis will be roughly perpendicular to the axes of the camera sensors. 

In the case of our pivot calibration, we would expect an increase in the stylus tracking error within one axis to translate to increased variation of the calculated stylus tips within the same axis. By observing that the stylus axis appoximately aligns with the z axis, we can limit our estimation of the probe's pose relative to the tracker to 2 dimensions (xy plane). We would suspect that the axis that is perpendicular to the stylus and shows the greatest variation in tip position corresponds with the axis perpendicular to the plane of the reflectors. Furthermore, since the location of the tracker relative to the pivot dimple did not change between our 1st and 2nd capture, we would expect the same axis to show a greater variation of points in both trials. This was nicely illustrated in both plots of our tip positions, where we found that the X positions of points tended to fall outside the 95% confidence interval more than in the other 2 dimensions (see additionalImage1 and additionalImage2 for reference). This axis is also perpendicular to the axis of the stylus, which matches our expectations. We could therefore expect the plane of the reflectors to lie roughly in the Y axis of the tracker, since this is both perpendicular to the axis of the probe (z-axis) and perpendicular to the axis of greatest variation (x-axis).

------------------------------------------------------
#6 - best fit vectors and standard deviation

Calibration_0:
best fit sphere radius: 188.3263
tip position in world CS: (129.02, -71.8028, -1420.23)
tip stdev in world CS:    (9.34618, 11.807, 16.4514)

Calibration_0_ransac:
% of inlier points: 90.74%
best fit sphere radius: 158.7604
tip position in world CS: (105.898, -67.1711, -1435.87)
tip stdev in world CS:    (0.797739, 1.20355, 3.22658)

Calibration_1:
best fit sphere radius: 160.0357
tip position in world CS: (121.005, -107.73, -1242.97)
tip stdev in world CS:    (0.387319, 0.518032, 0.69064)

Calibration_1_ransac:
% of inlier points: 100%
best fit sphere radius: 160.0357
tip position in world CS: (121.005, -107.73, -1242.97)
tip stdev in world CS:    (0.387319, 0.518032, 0.69064)

Calibration_2:
best fit sphere radius: 118.0384
tip position in world CS: (82.2916, -107.784, -1232.77)
tip stdev in world CS:    (11.3454, 19.5504, 15.6306)

Calibration_2_ransac:
% of inlier points: 98.38%
best fit sphere radius: 159.0171
tip position in world CS: (119.991, -107.898, -1242.69)
tip stdev in world CS:    (0.821695, 0.418305, 0.495829)