CISC 472 - Assignment 3 pt. 1


Students: Colton Barr 14cb49	Lauren Yates 15ly1

--------------------------------------------------

Answers and Results

2. ICP Results (without randomizing initialization)

knee1.csv
Rotation: [0.0006 0.0724 0.9974 ; -0.6632 0.7465 -0.0537 ; 0.7484 0.6615 -0.0484]
Translation: [1446.2866 -482.4646 -585.7567]
Number of Iterations: 15
RMSE: 19.29
Resulting Alignment: knee1_ICP_result.png

knee2.csv
Rotation: [0.5271 0.4025 0.7484 ; 0.4362 0.6277 -0.6448 ; -0.7293 0.6664 0.1553]
Translation: [1003.1852 -1328.4752 -276.3616]
Number of Iterations: 21
RMSE: 20.98
Resulting Alignment: knee2_ICP_result.png

4. ICP Results (with randomizing initialization)
knee1.csv
Rotation: [0.3000 0.9369 -0.1793; -0.8532 0.3476 0.3888; 0.4266 0.0364 0.9037]
Translation: [-479.9256 324.7711 992.4030]
Number of Iterations: 53
RMSE: 5.5336
Resulting Alignment: knee1_ICP_randomized_result.png

knee2.csv
Rotation: [-0.4270 -0.8908 0.1552; -0.8548 0.4536 0.2521; 0.2949 0.0250 0.9552]
Translation: [447.5441 105.4080 1076.7813]
Number of Iterations: 26
RMSE: 4.8988
Resulting Alignment: knee2_ICP_randomized_result.png

5. On average, the multi-attempt version of ICP worked better than the ICP without randomizing initialization. Rather than generate a truly random initial translation vector, we used our knowledge of the observed points' approximate locations on the knee to restrict the range of initial translations. This yielded a starting translation that was random but still reasonably close to the relevant part of the knee. Randomizing the initial rotation matrix also helped by generating unique alignments between the observed points and the knee, and in doing so the algorithm explored more of the search space and avoided local minima. The main benefit of multi-attempt ICP is the ability to better approximate an initial positioning for the points to the model, thus increasing the probability of having a more accurate alignment. However, the multi-attempt version of ICP is time consuming and typically requires several complete executions of the algorithm to reach a low RMSE solution. In addition, restricting the initial translations to a region of interest requires information about the target alignment (i.e. where the points are expected to lie on the model). 
