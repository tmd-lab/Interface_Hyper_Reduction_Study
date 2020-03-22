1

   Abaqus 3DEXPERIENCE R2019x                                  Date 19-Mar-2020   Time 00:25:31
   For use by WILLIAM MARSH RICE UNIVERSITY INC under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    1     INCREMENT     1     STEP TIME    0.00    
 


                        S T E P       1     M A T R I X   G E N E R A T E
 
  THE FOLLOWING GLOBAL MATRICES WILL BE GENERATED AND WRITTEN TO Job-3_X1.sim

    STIFFNESS
    MASS
 
  THE MASS MATRIX PROPERTIES: CONSISTENT MASS MATRIX
 
  IF MULTY-POINT CONSTRAINTS ARE SPECIFIED IN THE MODEL, THEY WILL BE APPLIED
 
  THE FOLLOWING GLOBAL MATRICES WILL BE WRITTEN TO THE FOLLOWING TEXT FILES

    STIFFNESS          TO Job-3_STIF1.mtx
    MASS               TO Job-3_MASS1.mtx
 
  MATRIX OUTPUT FORMAT - MATRIX INPUT
1

   Abaqus 3DEXPERIENCE R2019x                                  Date 19-Mar-2020   Time 00:25:32
   For use by WILLIAM MARSH RICE UNIVERSITY INC under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    2     INCREMENT     1     STEP TIME    0.00    


                        S T E P       2     C A L C U L A T I O N   O F   E I G E N V A L U E S

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

     CHARACTERISTIC ELEMENT LENGTH      2.00    

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     CHECK POINT   START OF SOLVER    
  
 COMPUTER PRECISION USED BY ABAQUS       2.22045e-16
 CONVERGENCE CRITERION FOR LANCZOS RUN   1.81899e-12
  
 NUMBER OF EQUATIONS                                                          18

 ***WARNING: ONLY 18 EIGENMODES WILL BE EXTRACTED BECAUSE THERE ARE NO MORE 
             UNRESTRAINED DEGREES OF FREEDOM
  
 NUMBER OF UNRESTRAINED DEGREES OF FREEDOM                                    18
 HEURISTIC SCALE FOR THE FIRST NONZERO EIGENVALUES (PROBLEM SCALE)   2.30200e+03
 INITIAL SHIFT VALUE                                                -2.30200e+03
  
 NUMBER OF CPUS USED BY DIRECT SOLVER     1
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 BLOCK SIZE WAS REDUCED. NEW BLOCK SIZE     1
  
 COMPUTATIONAL INTERVAL FOR THE INITIAL LANCZOS RUN
 
 EIGENVALUES NEEDED TO THE RIGHT OF SHIFT VALUE   18
 EIGENVALUES NEEDED TO THE LEFT OF SHIFT VALUE     0
 RIGHT BOUNDARY OF THE COMPUTATIONAL INTERVAL   3.95574e+37
 LEFT BOUNDARY OF THE COMPUTATIONAL INTERVAL   -3.95574e+37
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    9
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    1
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.52178e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -2.30200e+03
  
 SHIFTING TO THE RIGHT DIRECTION
 NEW SHIFT           4.52178e+07
  
 NUMBER OF NEGATIVE EIGENVALUES    18
 NUMBER OF POSITIVE EIGENVALUES     0
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: EXPAND THE TRUST INTERVAL TO THE RIGHT
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     12
 RIGHT BOUNDARY OF THE INTERVAL   4.52178e+07
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     12     0       |
 |  4.52178e+07  4.52178e+07 |  4.52178e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    2
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     12     0       |
 |  4.52178e+07  4.52178e+07 |  4.52178e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.52178e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           2.26078e+07
  
 NUMBER OF NEGATIVE EIGENVALUES    17
 NUMBER OF POSITIVE EIGENVALUES     1
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    17
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     12
 RIGHT BOUNDARY OF THE INTERVAL   2.26078e+07
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    2
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         2.26078e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.13027e+07
  
 NUMBER OF NEGATIVE EIGENVALUES    14
 NUMBER OF POSITIVE EIGENVALUES     4
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    14
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     1
 EIGENVALUES NEEDED TO THE LEFT     11
 RIGHT BOUNDARY OF THE INTERVAL   2.26078e+07
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.52178e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           5.65021e+06
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    12
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     1
 EIGENVALUES NEEDED TO THE LEFT     10
 RIGHT BOUNDARY OF THE INTERVAL   1.13027e+07
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  5.65021e+06  5.65021e+06 |  5.65021e+06  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  5.65021e+06  5.65021e+06 |  5.65021e+06  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.52178e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           2.82396e+06
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    12
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     10
 RIGHT BOUNDARY OF THE INTERVAL   2.82396e+06
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  2.82396e+06  5.65021e+06 |  2.82396e+06  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  2.82396e+06  5.65021e+06 |  2.82396e+06  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         2.82396e+06
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.41083e+06
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    12
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     10
 RIGHT BOUNDARY OF THE INTERVAL   1.41083e+06
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.41083e+06  5.65021e+06 |  1.41083e+06  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.41083e+06  5.65021e+06 |  1.41083e+06  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.41083e+06
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           7.04262e+05
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    12
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     10
 RIGHT BOUNDARY OF THE INTERVAL   7.04262e+05
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  7.04262e+05  5.65021e+06 |  7.04262e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  7.04262e+05  5.65021e+06 |  7.04262e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         7.04262e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           3.50980e+05
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    12
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     10
 RIGHT BOUNDARY OF THE INTERVAL   3.50980e+05
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  3.50980e+05  5.65021e+06 |  3.50980e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  3.50980e+05  5.65021e+06 |  3.50980e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         3.50980e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.74339e+05
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    12
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     10
 RIGHT BOUNDARY OF THE INTERVAL   1.74339e+05
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.74339e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           8.60185e+04
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     4
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   1.74339e+05
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  8.60185e+04  8.60185e+04 |  8.60185e+04  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  8.60185e+04  8.60185e+04 |  8.60185e+04  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.52178e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           4.18583e+04
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   4.18583e+04
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  4.18583e+04  8.60185e+04 |  4.18583e+04  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  4.18583e+04  8.60185e+04 |  4.18583e+04  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.18583e+04
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.97781e+04
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   1.97781e+04
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.97781e+04  8.60185e+04 |  1.97781e+04  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.97781e+04  8.60185e+04 |  1.97781e+04  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.97781e+04
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           8.73807e+03
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   8.73807e+03
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  8.73807e+03  8.60185e+04 |  8.73807e+03  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  8.73807e+03  8.60185e+04 |  8.73807e+03  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         8.73807e+03
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           3.21803e+03
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   3.21803e+03
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  3.21803e+03  8.60185e+04 |  3.21803e+03  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  3.21803e+03  8.60185e+04 |  3.21803e+03  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         3.21803e+03
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           4.58015e+02
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   4.58015e+02
 LEFT BOUNDARY OF THE INTERVAL   -2.30200e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  4.58015e+02  8.60185e+04 |  4.58015e+02  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     6 |      0    18       |
 | -3.95574e+37 -2.30200e+03 | -3.95574e+37 -2.30200e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  4.58015e+02  8.60185e+04 |  4.58015e+02  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.58015e+02
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -9.21993e+02
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   4.58015e+02
 LEFT BOUNDARY OF THE INTERVAL   -9.21993e+02
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  4.58015e+02  8.60185e+04 |  4.58015e+02  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -9.21993e+02 | -3.95574e+37 -9.21993e+02 | -3.95574e+37  4.58015e+02 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  4.58015e+02  8.60185e+04 |  4.58015e+02  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -9.21993e+02 | -3.95574e+37 -9.21993e+02 | -3.95574e+37  4.58015e+02 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.52178e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -9.21993e+02
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -2.31989e+02
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   4.58015e+02
 LEFT BOUNDARY OF THE INTERVAL   -2.31989e+02
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  4.58015e+02  8.60185e+04 |  4.58015e+02  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.31989e+02 | -3.95574e+37 -2.31989e+02 | -3.95574e+37  4.58015e+02 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  4.58015e+02  8.60185e+04 |  4.58015e+02  8.60185e+04 | -2.30200e+03  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.31989e+02 | -3.95574e+37 -2.31989e+02 | -3.95574e+37  4.58015e+02 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.52178e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -2.31989e+02
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.13013e+02
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   1.13013e+02
 LEFT BOUNDARY OF THE INTERVAL   -2.31989e+02
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.31989e+02 | -3.95574e+37 -2.31989e+02 | -3.95574e+37  4.58015e+02 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.13013e+02  8.60185e+04 |  1.13013e+02  8.60185e+04 | -2.31989e+02  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.31989e+02 | -3.95574e+37 -2.31989e+02 | -3.95574e+37  4.58015e+02 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.13013e+02  8.60185e+04 |  1.13013e+02  8.60185e+04 | -2.31989e+02  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.13013e+02
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -5.94879e+01
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.13013e+02
 LEFT BOUNDARY OF THE INTERVAL   -5.94879e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.13013e+02  8.60185e+04 |  1.13013e+02  8.60185e+04 | -2.31989e+02  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -5.94879e+01 | -3.95574e+37 -5.94879e+01 | -3.95574e+37  1.13013e+02 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.13013e+02  8.60185e+04 |  1.13013e+02  8.60185e+04 | -2.31989e+02  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -5.94879e+01 | -3.95574e+37 -5.94879e+01 | -3.95574e+37  1.13013e+02 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.52178e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -5.94879e+01
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           2.67626e+01
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   2.67626e+01
 LEFT BOUNDARY OF THE INTERVAL   -5.94879e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -5.94879e+01 | -3.95574e+37 -5.94879e+01 | -3.95574e+37  1.13013e+02 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  2.67626e+01  8.60185e+04 |  2.67626e+01  8.60185e+04 | -5.94879e+01  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -5.94879e+01 | -3.95574e+37 -5.94879e+01 | -3.95574e+37  1.13013e+02 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  2.67626e+01  8.60185e+04 |  2.67626e+01  8.60185e+04 | -5.94879e+01  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         2.67626e+01
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -1.63627e+01
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   2.67626e+01
 LEFT BOUNDARY OF THE INTERVAL   -1.63627e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  2.67626e+01  8.60185e+04 |  2.67626e+01  8.60185e+04 | -5.94879e+01  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.63627e+01 | -3.95574e+37 -1.63627e+01 | -3.95574e+37  2.67626e+01 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  2.67626e+01  8.60185e+04 |  2.67626e+01  8.60185e+04 | -5.94879e+01  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.63627e+01 | -3.95574e+37 -1.63627e+01 | -3.95574e+37  2.67626e+01 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.52178e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -1.63627e+01
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           5.19998e+00
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   5.19998e+00
 LEFT BOUNDARY OF THE INTERVAL   -1.63627e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.63627e+01 | -3.95574e+37 -1.63627e+01 | -3.95574e+37  2.67626e+01 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  5.19998e+00  8.60185e+04 |  5.19998e+00  8.60185e+04 | -1.63627e+01  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.63627e+01 | -3.95574e+37 -1.63627e+01 | -3.95574e+37  2.67626e+01 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  5.19998e+00  8.60185e+04 |  5.19998e+00  8.60185e+04 | -1.63627e+01  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         5.19998e+00
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -5.58134e+00
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   5.19998e+00
 LEFT BOUNDARY OF THE INTERVAL   -5.58134e+00
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  5.19998e+00  8.60185e+04 |  5.19998e+00  8.60185e+04 | -1.63627e+01  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -5.58134e+00 | -3.95574e+37 -5.58134e+00 | -3.95574e+37  5.19998e+00 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  5.19998e+00  8.60185e+04 |  5.19998e+00  8.60185e+04 | -1.63627e+01  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -5.58134e+00 | -3.95574e+37 -5.58134e+00 | -3.95574e+37  5.19998e+00 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.52178e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -5.58134e+00
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -1.90682e-01
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   5.19998e+00
 LEFT BOUNDARY OF THE INTERVAL   -1.90682e-01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  5.19998e+00  8.60185e+04 |  5.19998e+00  8.60185e+04 | -1.63627e+01  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.90682e-01 | -3.95574e+37 -1.90682e-01 | -3.95574e+37  5.19998e+00 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  5.19998e+00  8.60185e+04 |  5.19998e+00  8.60185e+04 | -1.63627e+01  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.90682e-01 | -3.95574e+37 -1.90682e-01 | -3.95574e+37  5.19998e+00 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.52178e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -1.90682e-01
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           2.50465e+00
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   2.50465e+00
 LEFT BOUNDARY OF THE INTERVAL   -1.90682e-01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.90682e-01 | -3.95574e+37 -1.90682e-01 | -3.95574e+37  5.19998e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  2.50465e+00  8.60185e+04 |  2.50465e+00  8.60185e+04 | -1.90682e-01  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.90682e-01 | -3.95574e+37 -1.90682e-01 | -3.95574e+37  5.19998e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  2.50465e+00  8.60185e+04 |  2.50465e+00  8.60185e+04 | -1.90682e-01  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         2.50465e+00
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.15698e+00
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   1.15698e+00
 LEFT BOUNDARY OF THE INTERVAL   -1.90682e-01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.90682e-01 | -3.95574e+37 -1.90682e-01 | -3.95574e+37  5.19998e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.15698e+00  8.60185e+04 |  1.15698e+00  8.60185e+04 | -1.90682e-01  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.90682e-01 | -3.95574e+37 -1.90682e-01 | -3.95574e+37  5.19998e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.15698e+00  8.60185e+04 |  1.15698e+00  8.60185e+04 | -1.90682e-01  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.15698e+00
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           4.83150e-01
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   4.83150e-01
 LEFT BOUNDARY OF THE INTERVAL   -1.90682e-01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.90682e-01 | -3.95574e+37 -1.90682e-01 | -3.95574e+37  5.19998e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  4.83150e-01  8.60185e+04 |  4.83150e-01  8.60185e+04 | -1.90682e-01  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.90682e-01 | -3.95574e+37 -1.90682e-01 | -3.95574e+37  5.19998e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  4.83150e-01  8.60185e+04 |  4.83150e-01  8.60185e+04 | -1.90682e-01  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.83150e-01
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.46234e-01
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   1.46234e-01
 LEFT BOUNDARY OF THE INTERVAL   -1.90682e-01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.90682e-01 | -3.95574e+37 -1.90682e-01 | -3.95574e+37  5.19998e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.46234e-01  8.60185e+04 |  1.46234e-01  8.60185e+04 | -1.90682e-01  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.90682e-01 | -3.95574e+37 -1.90682e-01 | -3.95574e+37  5.19998e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.46234e-01  8.60185e+04 |  1.46234e-01  8.60185e+04 | -1.90682e-01  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.46234e-01
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -2.22238e-02
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.46234e-01
 LEFT BOUNDARY OF THE INTERVAL   -2.22238e-02
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.46234e-01  8.60185e+04 |  1.46234e-01  8.60185e+04 | -1.90682e-01  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.22238e-02 | -3.95574e+37 -2.22238e-02 | -3.95574e+37  1.46234e-01 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   10
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.46234e-01  8.60185e+04 |  1.46234e-01  8.60185e+04 | -1.90682e-01  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.22238e-02 | -3.95574e+37 -2.22238e-02 | -3.95574e+37  1.46234e-01 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.52178e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -2.22238e-02
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           6.20052e-02
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   6.20052e-02
 LEFT BOUNDARY OF THE INTERVAL   -2.22238e-02
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.22238e-02 | -3.95574e+37 -2.22238e-02 | -3.95574e+37  1.46234e-01 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  6.20052e-02  8.60185e+04 |  6.20052e-02  8.60185e+04 | -2.22238e-02  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.22238e-02 | -3.95574e+37 -2.22238e-02 | -3.95574e+37  1.46234e-01 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  6.20052e-02  8.60185e+04 |  6.20052e-02  8.60185e+04 | -2.22238e-02  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.20052e-02
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.98907e-02
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   1.98907e-02
 LEFT BOUNDARY OF THE INTERVAL   -2.22238e-02
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.22238e-02 | -3.95574e+37 -2.22238e-02 | -3.95574e+37  1.46234e-01 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.98907e-02  8.60185e+04 |  1.98907e-02  8.60185e+04 | -2.22238e-02  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   10
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    6
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.22238e-02 | -3.95574e+37 -2.22238e-02 | -3.95574e+37  1.46234e-01 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.98907e-02  8.60185e+04 |  1.98907e-02  8.60185e+04 | -2.22238e-02  1.74339e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.98907e-02
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.51948e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -1.16654e-03
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.98907e-02
 LEFT BOUNDARY OF THE INTERVAL   -1.16654e-03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     4       |
 |  1.98907e-02  8.60185e+04 |  1.98907e-02  8.60185e+04 | -2.22238e-02  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.16654e-03 | -3.95574e+37 -1.16654e-03 | -3.95574e+37  1.98907e-02 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   10
 ALL THE EIGENVALUES REQUESTED ON THIS LANCZOS RUN ARE FOUND
  
 NUMBER OF ACCEPTED EIGENVALUES   12
  
 ACCEPTED EIGENVALUES
  -1.10569e-09 -3.28102e-10  1.12720e-10  3.18892e-10  4.08315e-10  7.48101e-10
   9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |     10     1       |
 |  1.74339e+05  5.65021e+06 |  1.74339e+05  5.65021e+06 | -2.30200e+03  1.13027e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     0 |      0     4       |
 | -3.95574e+37  8.60185e+04 | -3.95574e+37  8.60185e+04 | -3.95574e+37  1.74339e+05 |
 -------------------------------------------------------------------------------------
 SENTINEL WAS NOT INCLUDED INTO THE TRUST INTERVAL LIST. SHIFT  -1.16654e-03
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   12
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    6
 SHIFT RIGHT         4.58098e+03
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -1.16654e-03
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.22460e+05
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    12
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      4
 RIGHT BOUNDARY OF THE INTERVAL   1.22460e+05
 LEFT BOUNDARY OF THE INTERVAL    8.60185e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     0 |      0     4       |
 | -3.95574e+37  8.60185e+04 | -3.95574e+37  8.60185e+04 | -3.95574e+37  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      4     1       |
 |  1.22460e+05  5.65021e+06 |  1.22460e+05  5.65021e+06 |  8.60185e+04  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   13
  
 ACCEPTED EIGENVALUES
  -1.10569e-09 -3.28102e-10  1.12720e-10  3.18892e-10  4.08315e-10  7.48101e-10
   9.29783e+04  9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07  1.97239e+07
   4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     0 |      0     4       |
 | -3.95574e+37  8.60185e+04 | -3.95574e+37  8.60185e+04 | -3.95574e+37  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      4     1       |
 |  1.22460e+05  5.65021e+06 |  1.22460e+05  5.65021e+06 |  8.60185e+04  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   12
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.22460e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     1
 SHIFT LEFT          6.37917e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           1.22460e+05
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      3
 RIGHT BOUNDARY OF THE INTERVAL   1.22460e+05
 LEFT BOUNDARY OF THE INTERVAL    8.60185e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     0 |      0     4       |
 | -3.95574e+37  8.60185e+04 | -3.95574e+37  8.60185e+04 | -3.95574e+37  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     1     0 |      4     1       |
 |  1.22460e+05  5.65021e+06 |  1.22460e+05  5.65021e+06 |  8.60185e+04  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   14
  
 ACCEPTED EIGENVALUES
  -1.10569e-09 -3.28102e-10  1.12720e-10  3.18892e-10  4.08315e-10  7.48101e-10
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.90842e+06  1.40233e+07
   1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     0 |      0     4       |
 | -3.95574e+37  8.60185e+04 | -3.95574e+37  8.60185e+04 | -3.95574e+37  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     1     0 |      4     1       |
 |  1.22460e+05  5.65021e+06 |  1.22460e+05  5.65021e+06 |  8.60185e+04  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   12
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.22460e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     1
 SHIFT LEFT          6.37917e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           1.22460e+05
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      2
 RIGHT BOUNDARY OF THE INTERVAL   1.22460e+05
 LEFT BOUNDARY OF THE INTERVAL    8.60185e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     0 |      0     4       |
 | -3.95574e+37  8.60185e+04 | -3.95574e+37  8.60185e+04 | -3.95574e+37  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     1     0 |      4     1       |
 |  1.22460e+05  5.65021e+06 |  1.22460e+05  5.65021e+06 |  8.60185e+04  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   15
  
 ACCEPTED EIGENVALUES
  -1.10569e-09 -3.28102e-10  1.12720e-10  3.18892e-10  4.08315e-10  7.48101e-10
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.90842e+06
   1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     0 |      0     4       |
 | -3.95574e+37  8.60185e+04 | -3.95574e+37  8.60185e+04 | -3.95574e+37  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     1     0 |      4     1       |
 |  1.22460e+05  5.65021e+06 |  1.22460e+05  5.65021e+06 |  8.60185e+04  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   12
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.22460e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     1
 SHIFT LEFT          6.37917e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           1.22460e+05
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      1
 RIGHT BOUNDARY OF THE INTERVAL   1.22460e+05
 LEFT BOUNDARY OF THE INTERVAL    8.60185e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     0 |      0     4       |
 | -3.95574e+37  8.60185e+04 | -3.95574e+37  8.60185e+04 | -3.95574e+37  1.74339e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     1     0 |      4     1       |
 |  1.22460e+05  5.65021e+06 |  1.22460e+05  5.65021e+06 |  8.60185e+04  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES REQUESTED ON THIS LANCZOS RUN ARE FOUND
  
 NUMBER OF ACCEPTED EIGENVALUES   16
  
 ACCEPTED EIGENVALUES
  -1.10569e-09 -3.28102e-10  1.12720e-10  3.18892e-10  4.08315e-10  7.48101e-10
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  +  |     12 |        0       12 |     0     0 |     0     0 |      0     1       |
 | -3.95574e+37  5.65021e+06 | -3.95574e+37  5.65021e+06 | -3.95574e+37  1.13027e+07 |
 -------------------------------------------------------------------------------------
 SENTINEL WAS NOT INCLUDED INTO THE TRUST INTERVAL LIST. SHIFT   1.22460e+05
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL    6
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.22460e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     1
 SHIFT LEFT          6.37917e+04
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           7.99142e+06
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    12
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     1
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.13027e+07
 LEFT BOUNDARY OF THE INTERVAL    7.99142e+06
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       14       14 |     0     0 |     0     0 |     11     1       |
 |  1.13027e+07  1.13027e+07 |  1.13027e+07  1.13027e+07 | -2.30200e+03  2.26078e+07 |
 -------------------------------------------------------------------------------------
 |  +  |     12 |        0       12 |     0     0 |     0     0 |      0     1       |
 | -3.95574e+37  7.99142e+06 | -3.95574e+37  7.99142e+06 | -3.95574e+37  1.13027e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    6
 ALL THE EIGENVALUES REQUESTED ON THIS LANCZOS RUN ARE FOUND
  
 NUMBER OF ACCEPTED EIGENVALUES   17
  
 ACCEPTED EIGENVALUES
  -1.10569e-09 -3.28102e-10  1.12720e-10  3.18892e-10  4.08315e-10  7.48101e-10
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.90842e+06  9.90842e+06  1.40233e+07  1.97239e+07  4.51948e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    2
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  +  |     14 |        0       14 |     0     0 |     0     0 |      0     1       |
 | -3.95574e+37  1.13027e+07 | -3.95574e+37  1.13027e+07 | -3.95574e+37  2.26078e+07 |
 -------------------------------------------------------------------------------------
 SENTINEL WAS NOT INCLUDED INTO THE TRUST INTERVAL LIST. SHIFT   7.99142e+06
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL    4
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.18063e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.99142e+06
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.59853e+07
  
 NUMBER OF NEGATIVE EIGENVALUES    15
 NUMBER OF POSITIVE EIGENVALUES     3
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    15
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     1
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   2.26078e+07
 LEFT BOUNDARY OF THE INTERVAL    1.59853e+07
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       17       18 |     0     0 |     0     0 |     12     0       |
 |  2.26078e+07  4.52178e+07 |  2.26078e+07  4.52178e+07 | -2.30200e+03  4.52178e+07 |
 -------------------------------------------------------------------------------------
 |  +  |     15 |        0       15 |     0     0 |     0     0 |      0     1       |
 | -3.95574e+37  1.59853e+07 | -3.95574e+37  1.59853e+07 | -3.95574e+37  2.26078e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    6
 ALL THE EIGENVALUES REQUESTED ON THIS LANCZOS RUN ARE FOUND
  
 NUMBER OF ACCEPTED EIGENVALUES   18
  
 ACCEPTED EIGENVALUES
  -1.10569e-09 -3.28102e-10  1.12720e-10  3.18892e-10  4.08315e-10  7.48101e-10
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.90842e+06  9.90842e+06  1.40233e+07  1.97239e+07  1.97239e+07  4.51948e+07
  
 THE REQUIRED NUMBER OF EIGENVALUES IS FOUND
  
 FINAL TRUST INTERVALS AFTER SUCCESSFUL LANCZOS RUN
  
 NUMBER OF TRUST INTERVALS    1
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |     18 |        0       18 |     0     0 |     0     0 |      0     0       |
 | -3.95574e+37  4.52178e+07 | -3.95574e+37  4.52178e+07 | -3.95574e+37  4.52178e+07 |
 -------------------------------------------------------------------------------------
  
 FIRST EIGENVALUE REQUESTED BY THE USER     1
  
 EIGENVALUES REQUESTED BY THE USER
  -1.10569e-09 -3.28102e-10  1.12720e-10  3.18892e-10  4.08315e-10  7.48101e-10
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.90842e+06  9.90842e+06  1.40233e+07  1.97239e+07  1.97239e+07  4.51948e+07
  
 EIGENVALUE EXTRACTION USING THE LANCZOS EIGENSOLVER IS COMPLETE
 NUMBER OF LANCZOS RUNS       :  39
 NUMBER OF LANCZOS STEPS      :  92
 NUMBER OF FACTORIZATIONS     :  36
 NUMBER OF BACK SUBSTITUTIONS :  99

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.27000    
       ELAPSED SYSTEM TIME (SEC)    =  0.10000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.37000    
       ELAPSED WALLCLOCK TIME (SEC) =          0


          THE ANALYSIS HAS BEEN COMPLETED



     ANALYSIS SUMMARY:
     TOTAL OF          2  INCREMENTS
                       0  CUTBACKS IN AUTOMATIC INCREMENTATION
                       0  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                     135  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                      36  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       0  WARNING MESSAGES DURING USER INPUT PROCESSING
                       1  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                       0  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       0  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =   1.0000    
       SYSTEM TIME (SEC)    =  0.13000    
       TOTAL CPU TIME (SEC) =   1.1300    
       WALLCLOCK TIME (SEC) =          2
