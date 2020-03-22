1

   Abaqus 3DEXPERIENCE R2019x                                  Date 17-Mar-2020   Time 22:55:36
   For use by WILLIAM MARSH RICE UNIVERSITY INC under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    1     INCREMENT     1     STEP TIME    0.00    


                        S T E P       1     C A L C U L A T I O N   O F   E I G E N V A L U E S

                                            F O R   N A T U R A L   F R E Q U E N C I E S



                                                                                          

     THE LANCZOS EIGENSOLVER IS USED FOR THIS ANALYSIS
     SYSTEM DATA WILL BE WRITTEN TO THE .SIM FILE

     Abaqus WILL COMPUTE UNCOUPLED
      STRUCTURAL AND ACOUSTIC MODES
     NUMBER OF EIGENVALUES                     60
     HIGHEST FREQUENCY OF INTEREST                  1.00000E+18
     MAXIMUM NUMBER OF STEPS WITHIN RUN               35
     BLOCK SIZE FOR LANCZOS PROCEDURE                  7
     THE EIGENVECTORS ARE SCALED SO THAT 
      THE GENERALIZED MASS IN EACH VECTOR IS UNITY

     THIS IS A LINEAR PERTURBATION STEP.
     ALL LOADS ARE DEFINED AS CHANGE IN LOAD TO THE REFERENCE STATE

     EXTRAPOLATION WILL NOT BE USED

     CHARACTERISTIC ELEMENT LENGTH     1.906E-03

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     CHECK POINT   START OF SOLVER    
  
 COMPUTER PRECISION USED BY ABAQUS       2.22045e-16
 CONVERGENCE CRITERION FOR LANCZOS RUN   1.81899e-12
  
 NUMBER OF EQUATIONS                                                        6180
 NUMBER OF UNRESTRAINED DEGREES OF FREEDOM                                  6120
 HEURISTIC SCALE FOR THE FIRST NONZERO EIGENVALUES (PROBLEM SCALE)   2.43957e+05
 INITIAL SHIFT VALUE                                                 0.00000e+00
  
 NUMBER OF CPUS USED BY DIRECT SOLVER     1
  
 NUMBER OF NEGATIVE EIGENVALUES    35
 NUMBER OF POSITIVE EIGENVALUES  6145
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS    48
  
 COMPUTATIONAL INTERVAL FOR THE INITIAL LANCZOS RUN
 
 EIGENVALUES NEEDED TO THE RIGHT OF SHIFT VALUE   60
 EIGENVALUES NEEDED TO THE LEFT OF SHIFT VALUE     0
 RIGHT BOUNDARY OF THE COMPUTATIONAL INTERVAL   3.95574e+37
 LEFT BOUNDARY OF THE COMPUTATIONAL INTERVAL   -3.95574e+37
 INITIAL BLOCK CONTAINS LINEARLY DEPENDENT VECTORS. REDUCING THE BLOCK SIZE ACCORDINGLY
  
 NEW BLOCK SIZE               6
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
  -3.51625e-04 -1.98793e-04 -1.16153e-04 -5.15998e-05 -4.90841e-05  1.16133e-04
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    1
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |      5 |       30       35 |     0     0 |     0     0 |      0    60       |
 | -3.95574e+37  0.00000e+00 | -3.95574e+37  0.00000e+00 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   55
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         4.85474e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          0.00000e+00
  
 SHIFTING TO THE RIGHT DIRECTION
 NEW SHIFT           4.85474e+05
  
 NUMBER OF NEGATIVE EIGENVALUES    36
 NUMBER OF POSITIVE EIGENVALUES  6144
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS    48
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: EXPAND THE TRUST INTERVAL TO THE RIGHT
 EIGENVALUES NEEDED TO THE RIGHT    54
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   3.95574e+37
 LEFT BOUNDARY OF THE INTERVAL    4.85474e+05
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |      6 |       30       36 |     0     0 |     0     0 |      0    54       |
 | -3.95574e+37  4.85474e+05 | -3.95574e+37  4.85474e+05 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   15
 TERMINATION BASED ON THE COST-PER-EIGENVALUE CRITERION
  
 NUMBER OF ACCEPTED EIGENVALUES   27
  
 ACCEPTED EIGENVALUES
  -3.51625e-04 -1.98793e-04 -1.16153e-04 -5.15998e-05 -4.90841e-05  1.16133e-04
   2.23970e+07  1.73515e+08  1.76275e+08  6.67313e+08  7.48609e+08  1.84684e+09
   1.93491e+09  2.46231e+09  3.89610e+09  4.23029e+09  6.48528e+09  8.11004e+09
   1.03463e+10  1.07313e+10  1.15431e+10  1.25278e+10  1.31159e+10  1.45108e+10
   1.52435e+10  1.68098e+10  1.74996e+10
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    1
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |      6 |       30       36 |     0     0 |     0     0 |      0    54       |
 | -3.95574e+37  4.85474e+05 | -3.95574e+37  4.85474e+05 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   54
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT   21
 SHIFT RIGHT         6.43379e+10
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          4.85474e+05
  
 SHIFTING TO THE RIGHT DIRECTION
 NEW SHIFT           6.43379e+10
  
 NUMBER OF NEGATIVE EIGENVALUES   137
 NUMBER OF POSITIVE EIGENVALUES  6043
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    89
 NUMBER OF LAGRANGE MULTIPLIERS    48
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: EXPAND THE TRUST INTERVAL TO THE RIGHT
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     80
 RIGHT BOUNDARY OF THE INTERVAL   6.43379e+10
 LEFT BOUNDARY OF THE INTERVAL    1.74996e+10
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |      6 |       30       36 |     0    21 |     0    21 |      0    54       |
 | -3.95574e+37  4.85474e+05 | -3.95574e+37  1.74996e+10 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      137      137 |     0     0 |     0     0 |     80     0       |
 |  6.43379e+10  6.43379e+10 |  6.43379e+10  6.43379e+10 |  1.74996e+10  6.43379e+10 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   15
 TERMINATION BASED ON THE COST-PER-EIGENVALUE CRITERION
  
 NUMBER OF ACCEPTED EIGENVALUES   40
  
 ACCEPTED EIGENVALUES
  -3.51625e-04 -1.98793e-04 -1.16153e-04 -5.15998e-05 -4.90841e-05  1.16133e-04
   2.23970e+07  1.73515e+08  1.76275e+08  6.67313e+08  7.48609e+08  1.84684e+09
   1.93491e+09  2.46231e+09  3.89610e+09  4.23029e+09  6.48528e+09  8.11004e+09
   1.03463e+10  1.07313e+10  1.15431e+10  1.25278e+10  1.31159e+10  1.45108e+10
   1.52435e+10  1.68098e+10  1.74996e+10  6.29256e+10  6.30202e+10  6.31466e+10
   6.31538e+10  6.33487e+10  6.34254e+10  6.35439e+10  6.39571e+10  6.40551e+10
   6.40793e+10  6.40989e+10  6.41404e+10  6.41650e+10
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    2
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |      6 |       30       36 |     0    21 |     0    21 |      0    54       |
 | -3.95574e+37  4.85474e+05 | -3.95574e+37  1.74996e+10 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      137      137 |     0     0 |     0     0 |     80     0       |
 |  6.43379e+10  6.43379e+10 |  6.43379e+10  6.43379e+10 |  1.74996e+10  6.43379e+10 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   54
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.43379e+10
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT    13
 SHIFT LEFT          2.95146e+10
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           3.31839e+10
  
 NUMBER OF NEGATIVE EIGENVALUES    61
 NUMBER OF POSITIVE EIGENVALUES  6119
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    13
 NUMBER OF LAGRANGE MULTIPLIERS    48
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT    63
 EIGENVALUES NEEDED TO THE LEFT      4
 RIGHT BOUNDARY OF THE INTERVAL   6.29256e+10
 LEFT BOUNDARY OF THE INTERVAL    1.74996e+10
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |      6 |       30       36 |     0    21 |     0    21 |      0    54       |
 | -3.95574e+37  4.85474e+05 | -3.95574e+37  1.74996e+10 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      137      137 |    13     0 |    13     0 |     80     0       |
 |  6.43379e+10  6.43379e+10 |  6.29256e+10  6.43379e+10 |  1.74996e+10  6.43379e+10 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       61       61 |     0     0 |     0     0 |      4    63       |
 |  3.31839e+10  3.31839e+10 |  3.31839e+10  3.31839e+10 |  1.74996e+10  6.29256e+10 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   15
 TERMINATION BASED ON THE COST-PER-EIGENVALUE CRITERION
  
 NUMBER OF ACCEPTED EIGENVALUES   51
  
 ACCEPTED EIGENVALUES
  -3.51625e-04 -1.98793e-04 -1.16153e-04 -5.15998e-05 -4.90841e-05  1.16133e-04
   2.23970e+07  1.73515e+08  1.76275e+08  6.67313e+08  7.48609e+08  1.84684e+09
   1.93491e+09  2.46231e+09  3.89610e+09  4.23029e+09  6.48528e+09  8.11004e+09
   1.03463e+10  1.07313e+10  1.15431e+10  1.25278e+10  1.31159e+10  1.45108e+10
   1.52435e+10  1.68098e+10  1.74996e+10  2.29941e+10  2.32817e+10  2.53136e+10
   3.04742e+10  3.57552e+10  3.59679e+10  3.63658e+10  3.97545e+10  4.90667e+10
   4.97817e+10  5.09696e+10  6.29256e+10  6.30202e+10  6.31466e+10  6.31538e+10
   6.33487e+10  6.34254e+10  6.35439e+10  6.39571e+10  6.40551e+10  6.40793e+10
   6.40989e+10  6.41404e+10  6.41650e+10
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    2
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      137      137 |    13     0 |    13     0 |     80     0       |
 |  6.43379e+10  6.43379e+10 |  6.29256e+10  6.43379e+10 |  1.74996e+10  6.43379e+10 |
 -------------------------------------------------------------------------------------
 |  +  |     31 |       30       61 |     0     0 |     0     0 |      0    63       |
 | -3.95574e+37  3.31839e+10 | -3.95574e+37  3.31839e+10 | -3.95574e+37  6.29256e+10 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   29
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    7
 SHIFT RIGHT         6.85774e+10
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     4
 SHIFT LEFT         -1.63935e+09
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           5.69476e+10
  
 NUMBER OF NEGATIVE EIGENVALUES    70
 NUMBER OF POSITIVE EIGENVALUES  6110
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    22
 NUMBER OF LAGRANGE MULTIPLIERS    48
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT    54
 EIGENVALUES NEEDED TO THE LEFT      2
 RIGHT BOUNDARY OF THE INTERVAL   6.29256e+10
 LEFT BOUNDARY OF THE INTERVAL    5.09696e+10
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      137      137 |    13     0 |    13     0 |     80     0       |
 |  6.43379e+10  6.43379e+10 |  6.29256e+10  6.43379e+10 |  1.74996e+10  6.43379e+10 |
 -------------------------------------------------------------------------------------
 |  +  |     31 |       30       61 |     0     7 |     0     7 |      0    63       |
 | -3.95574e+37  3.31839e+10 | -3.95574e+37  5.09696e+10 | -3.95574e+37  6.29256e+10 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       70       70 |     0     0 |     0     0 |      2    54       |
 |  5.69476e+10  5.69476e+10 |  5.69476e+10  5.69476e+10 |  5.09696e+10  6.29256e+10 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   18
 TERMINATION BASED ON THE COST-PER-EIGENVALUE CRITERION
  
 NUMBER OF ACCEPTED EIGENVALUES   71
  
 ACCEPTED EIGENVALUES
  -3.51625e-04 -1.98793e-04 -1.16153e-04 -5.15998e-05 -4.90841e-05  1.16133e-04
   2.23970e+07  1.73515e+08  1.76275e+08  6.67313e+08  7.48609e+08  1.84684e+09
   1.93491e+09  2.46231e+09  3.89610e+09  4.23029e+09  6.48528e+09  8.11004e+09
   1.03463e+10  1.07313e+10  1.15431e+10  1.25278e+10  1.31159e+10  1.45108e+10
   1.52435e+10  1.68098e+10  1.74996e+10  2.29941e+10  2.32817e+10  2.53136e+10
   3.04742e+10  3.57552e+10  3.59679e+10  3.63658e+10  3.97545e+10  4.90667e+10
   4.97817e+10  5.09696e+10  5.45648e+10  5.54249e+10  5.77242e+10  5.77254e+10
   5.78185e+10  5.78223e+10  5.83368e+10  5.83401e+10  5.83415e+10  5.83431e+10
   5.83669e+10  5.83850e+10  5.92130e+10  5.92609e+10  5.97743e+10  5.98461e+10
   5.98953e+10  6.06560e+10  6.18601e+10  6.21524e+10  6.29256e+10  6.30202e+10
   6.31466e+10  6.31538e+10  6.33487e+10  6.34254e+10  6.35439e+10  6.39571e+10
   6.40551e+10  6.40793e+10  6.40989e+10  6.41404e+10  6.41650e+10
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    2
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      137      137 |    13     0 |    13     0 |     80     0       |
 |  6.43379e+10  6.43379e+10 |  6.29256e+10  6.43379e+10 |  1.74996e+10  6.43379e+10 |
 -------------------------------------------------------------------------------------
 |  +  |     40 |       30       70 |     0     0 |     0     0 |      0    54       |
 | -3.95574e+37  5.69476e+10 | -3.95574e+37  5.69476e+10 | -3.95574e+37  6.29256e+10 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   20
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         2.39156e+11
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     2
 SHIFT LEFT          2.21243e+10
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           5.99366e+10
  
 NUMBER OF NEGATIVE EIGENVALUES   112
 NUMBER OF POSITIVE EIGENVALUES  6068
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    64
 NUMBER OF LAGRANGE MULTIPLIERS    48
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     9
 EIGENVALUES NEEDED TO THE LEFT     27
 RIGHT BOUNDARY OF THE INTERVAL   6.29256e+10
 LEFT BOUNDARY OF THE INTERVAL    5.69476e+10
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      137      137 |    13     0 |    13     0 |     80     0       |
 |  6.43379e+10  6.43379e+10 |  6.29256e+10  6.43379e+10 |  1.74996e+10  6.43379e+10 |
 -------------------------------------------------------------------------------------
 |  +  |     40 |       30       70 |     0     0 |     0    18 |      0    54       |
 | -3.95574e+37  5.69476e+10 | -3.95574e+37  5.69476e+10 | -3.95574e+37  6.29256e+10 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      112      112 |     0     0 |     0     0 |     27     9       |
 |  5.99366e+10  5.99366e+10 |  5.99366e+10  5.99366e+10 |  5.69476e+10  6.29256e+10 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   15
 TERMINATION BASED ON THE COST-PER-EIGENVALUE CRITERION
  
 NUMBER OF ACCEPTED EIGENVALUES  101
  
 ACCEPTED EIGENVALUES
  -3.51625e-04 -1.98793e-04 -1.16153e-04 -5.15998e-05 -4.90841e-05  1.16133e-04
   2.23970e+07  1.73515e+08  1.76275e+08  6.67313e+08  7.48609e+08  1.84684e+09
   1.93491e+09  2.46231e+09  3.89610e+09  4.23029e+09  6.48528e+09  8.11004e+09
   1.03463e+10  1.07313e+10  1.15431e+10  1.25278e+10  1.31159e+10  1.45108e+10
   1.52435e+10  1.68098e+10  1.74996e+10  2.29941e+10  2.32817e+10  2.53136e+10
   3.04742e+10  3.57552e+10  3.59679e+10  3.63658e+10  3.97545e+10  4.90667e+10
   4.97817e+10  5.09696e+10  5.45648e+10  5.54249e+10  5.76467e+10  5.76467e+10
   5.76467e+10  5.76467e+10  5.76467e+10  5.76467e+10  5.77242e+10  5.77254e+10
   5.78185e+10  5.78223e+10  5.83368e+10  5.83401e+10  5.83415e+10  5.83431e+10
   5.83669e+10  5.83850e+10  5.92130e+10  5.92609e+10  5.95059e+10  5.95060e+10
   5.95061e+10  5.95061e+10  5.95061e+10  5.95061e+10  5.95061e+10  5.95061e+10
   5.95061e+10  5.95061e+10  5.95062e+10  5.95062e+10  5.97273e+10  5.97356e+10
   5.97743e+10  5.98461e+10  5.98953e+10  5.99193e+10  6.00581e+10  6.00747e+10
   6.02566e+10  6.02967e+10  6.03732e+10  6.04571e+10  6.05235e+10  6.06560e+10
   6.15910e+10  6.16128e+10  6.18601e+10  6.21524e+10  6.29256e+10  6.30202e+10
   6.31466e+10  6.31538e+10  6.33487e+10  6.34254e+10  6.35439e+10  6.39571e+10
   6.40551e+10  6.40793e+10  6.40989e+10  6.41404e+10  6.41650e+10
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    2
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |     40 |       30       70 |     0     0 |     0    18 |      0    54       |
 | -3.95574e+37  5.69476e+10 | -3.95574e+37  5.69476e+10 | -3.95574e+37  6.29256e+10 |
 -------------------------------------------------------------------------------------
 |  -  |     25 |      112      137 |     0     0 |     0     0 |     27     0       |
 |  5.99366e+10  6.43379e+10 |  5.99366e+10  6.43379e+10 |  5.69476e+10  6.43379e+10 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   20
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    9
 SHIFT RIGHT         9.47599e+10
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT    21
 SHIFT LEFT          2.51133e+10
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           5.99366e+10
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   5.99366e+10
 LEFT BOUNDARY OF THE INTERVAL    5.69476e+10
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |     40 |       30       70 |     0     0 |     0    18 |      0    54       |
 | -3.95574e+37  5.69476e+10 | -3.95574e+37  5.69476e+10 | -3.95574e+37  6.29256e+10 |
 -------------------------------------------------------------------------------------
 |  -  |     25 |      112      137 |     0     0 |    21     0 |     27     0       |
 |  5.99366e+10  6.43379e+10 |  5.99366e+10  6.43379e+10 |  5.69476e+10  6.43379e+10 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   17
 ALL THE EIGENVALUES REQUESTED ON THIS LANCZOS RUN ARE FOUND
  
 NUMBER OF ACCEPTED EIGENVALUES  107
  
 ACCEPTED EIGENVALUES
  -3.51625e-04 -1.98793e-04 -1.16153e-04 -5.15998e-05 -4.90841e-05  1.16133e-04
   2.23970e+07  1.73515e+08  1.76275e+08  6.67313e+08  7.48609e+08  1.84684e+09
   1.93491e+09  2.46231e+09  3.89610e+09  4.23029e+09  6.48528e+09  8.11004e+09
   1.03463e+10  1.07313e+10  1.15431e+10  1.25278e+10  1.31159e+10  1.45108e+10
   1.52435e+10  1.68098e+10  1.74996e+10  2.29941e+10  2.32817e+10  2.53136e+10
   3.04742e+10  3.57552e+10  3.59679e+10  3.63658e+10  3.97545e+10  4.90667e+10
   4.97817e+10  5.09696e+10  5.45648e+10  5.54249e+10  5.76467e+10  5.76467e+10
   5.76467e+10  5.76467e+10  5.76467e+10  5.76467e+10  5.76467e+10  5.76467e+10
   5.76467e+10  5.76467e+10  5.76467e+10  5.76467e+10  5.77242e+10  5.77254e+10
   5.78185e+10  5.78223e+10  5.83368e+10  5.83401e+10  5.83415e+10  5.83431e+10
   5.83669e+10  5.83850e+10  5.92130e+10  5.92609e+10  5.95059e+10  5.95060e+10
   5.95061e+10  5.95061e+10  5.95061e+10  5.95061e+10  5.95061e+10  5.95061e+10
   5.95061e+10  5.95061e+10  5.95062e+10  5.95062e+10  5.97273e+10  5.97356e+10
   5.97743e+10  5.98461e+10  5.98953e+10  5.99193e+10  6.00581e+10  6.00747e+10
   6.02566e+10  6.02967e+10  6.03732e+10  6.04571e+10  6.05235e+10  6.06560e+10
   6.15910e+10  6.16128e+10  6.18601e+10  6.21524e+10  6.29256e+10  6.30202e+10
   6.31466e+10  6.31538e+10  6.33487e+10  6.34254e+10  6.35439e+10  6.39571e+10
   6.40551e+10  6.40793e+10  6.40989e+10  6.41404e+10  6.41650e+10
  
 THE REQUIRED NUMBER OF EIGENVALUES IS FOUND
  
 FINAL TRUST INTERVALS AFTER SUCCESSFUL LANCZOS RUN
  
 NUMBER OF TRUST INTERVALS    1
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |    107 |       30      137 |     0     0 |     0     0 |      0     0       |
 | -3.95574e+37  6.43379e+10 | -3.95574e+37  6.43379e+10 | -3.95574e+37  6.43379e+10 |
 -------------------------------------------------------------------------------------
  
 FIRST EIGENVALUE REQUESTED BY THE USER     1
  
 EIGENVALUES REQUESTED BY THE USER
  -3.51625e-04 -1.98793e-04 -1.16153e-04 -5.15998e-05 -4.90841e-05  1.16133e-04
   2.23970e+07  1.73515e+08  1.76275e+08  6.67313e+08  7.48609e+08  1.84684e+09
   1.93491e+09  2.46231e+09  3.89610e+09  4.23029e+09  6.48528e+09  8.11004e+09
   1.03463e+10  1.07313e+10  1.15431e+10  1.25278e+10  1.31159e+10  1.45108e+10
   1.52435e+10  1.68098e+10  1.74996e+10  2.29941e+10  2.32817e+10  2.53136e+10
   3.04742e+10  3.57552e+10  3.59679e+10  3.63658e+10  3.97545e+10  4.90667e+10
   4.97817e+10  5.09696e+10  5.45648e+10  5.54249e+10  5.76467e+10  5.76467e+10
   5.76467e+10  5.76467e+10  5.76467e+10  5.76467e+10  5.76467e+10  5.76467e+10
   5.76467e+10  5.76467e+10  5.76467e+10  5.76467e+10  5.77242e+10  5.77254e+10
   5.78185e+10  5.78223e+10  5.83368e+10  5.83401e+10  5.83415e+10  5.83431e+10
  
 EIGENVALUE EXTRACTION USING THE LANCZOS EIGENSOLVER IS COMPLETE
 NUMBER OF LANCZOS RUNS       :   7
 NUMBER OF LANCZOS STEPS      :  96
 NUMBER OF FACTORIZATIONS     :   6
 NUMBER OF BACK SUBSTITUTIONS : 100

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   2.5100    
       ELAPSED SYSTEM TIME (SEC)    =  0.44000    
       ELAPSED TOTAL CPU TIME (SEC) =   2.9500    
       ELAPSED WALLCLOCK TIME (SEC) =          3


          THE ANALYSIS HAS BEEN COMPLETED



     ANALYSIS SUMMARY:
     TOTAL OF          1  INCREMENTS
                       0  CUTBACKS IN AUTOMATIC INCREMENTATION
                       0  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                     106  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                       6  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       5  WARNING MESSAGES DURING USER INPUT PROCESSING
                       0  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                       0  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       0  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =   3.1100    
       SYSTEM TIME (SEC)    =  0.48000    
       TOTAL CPU TIME (SEC) =   3.5900    
       WALLCLOCK TIME (SEC) =          4
