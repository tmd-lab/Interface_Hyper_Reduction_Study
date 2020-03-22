1

   Abaqus 3DEXPERIENCE R2019x                                  Date 19-Mar-2020   Time 00:19:20
   For use by WILLIAM MARSH RICE UNIVERSITY INC under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    1     INCREMENT     1     STEP TIME    0.00    
 


                        S T E P       1     M A T R I X   G E N E R A T E
 
  THE FOLLOWING GLOBAL MATRICES WILL BE GENERATED AND WRITTEN TO Job-2_X1.sim

    STIFFNESS
    MASS
 
  THE MASS MATRIX PROPERTIES: CONSISTENT MASS MATRIX
 
  IF MULTY-POINT CONSTRAINTS ARE SPECIFIED IN THE MODEL, THEY WILL BE APPLIED
 
  THE FOLLOWING GLOBAL MATRICES WILL BE WRITTEN TO THE FOLLOWING TEXT FILES

    STIFFNESS          TO Job-2_STIF1.mtx
    MASS               TO Job-2_MASS1.mtx
 
  MATRIX OUTPUT FORMAT - MATRIX INPUT
1

   Abaqus 3DEXPERIENCE R2019x                                  Date 19-Mar-2020   Time 00:19:22
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
  
 NUMBER OF EQUATIONS                                                          36

 ***WARNING: ONLY 24 EIGENMODES WILL BE EXTRACTED BECAUSE THERE ARE NO MORE 
             UNRESTRAINED DEGREES OF FREEDOM
  
 NUMBER OF UNRESTRAINED DEGREES OF FREEDOM                                    24
 HEURISTIC SCALE FOR THE FIRST NONZERO EIGENVALUES (PROBLEM SCALE)   9.48134e+03
 INITIAL SHIFT VALUE                                                -9.48134e+03
  
 NUMBER OF CPUS USED BY DIRECT SOLVER     1
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    30
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 BLOCK SIZE WAS REDUCED. NEW BLOCK SIZE     1
  
 COMPUTATIONAL INTERVAL FOR THE INITIAL LANCZOS RUN
 
 EIGENVALUES NEEDED TO THE RIGHT OF SHIFT VALUE   24
 EIGENVALUES NEEDED TO THE LEFT OF SHIFT VALUE     0
 RIGHT BOUNDARY OF THE COMPUTATIONAL INTERVAL   3.95574e+37
 LEFT BOUNDARY OF THE COMPUTATIONAL INTERVAL   -3.95574e+37
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   13
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    1
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.41974e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -9.48134e+03
  
 SHIFTING TO THE RIGHT DIRECTION
 NEW SHIFT           6.41974e+07
  
 NUMBER OF NEGATIVE EIGENVALUES    30
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    24
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: EXPAND THE TRUST INTERVAL TO THE RIGHT
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     16
 RIGHT BOUNDARY OF THE INTERVAL   6.41974e+07
 LEFT BOUNDARY OF THE INTERVAL   -9.48134e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       30       30 |     0     0 |     0     0 |     16     0       |
 |  6.41974e+07  6.41974e+07 |  6.41974e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    2
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       30       30 |     0     0 |     0     0 |     16     0       |
 |  6.41974e+07  6.41974e+07 |  6.41974e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.41974e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           3.20939e+07
  
 NUMBER OF NEGATIVE EIGENVALUES    29
 NUMBER OF POSITIVE EIGENVALUES     7
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    23
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     16
 RIGHT BOUNDARY OF THE INTERVAL   3.20939e+07
 LEFT BOUNDARY OF THE INTERVAL   -9.48134e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    2
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         3.20939e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.60422e+07
  
 NUMBER OF NEGATIVE EIGENVALUES    24
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     2
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   3.20939e+07
 LEFT BOUNDARY OF THE INTERVAL   -9.48134e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.60422e+07  1.60422e+07 |  1.60422e+07  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.60422e+07  1.60422e+07 |  1.60422e+07  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.41974e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           8.01638e+06
  
 NUMBER OF NEGATIVE EIGENVALUES    24
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   8.01638e+06
 LEFT BOUNDARY OF THE INTERVAL   -9.48134e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  8.01638e+06  1.60422e+07 |  8.01638e+06  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  8.01638e+06  1.60422e+07 |  8.01638e+06  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         8.01638e+06
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           4.00345e+06
  
 NUMBER OF NEGATIVE EIGENVALUES    24
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   4.00345e+06
 LEFT BOUNDARY OF THE INTERVAL   -9.48134e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  4.00345e+06  1.60422e+07 |  4.00345e+06  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  4.00345e+06  1.60422e+07 |  4.00345e+06  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.00345e+06
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.99698e+06
  
 NUMBER OF NEGATIVE EIGENVALUES    24
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   1.99698e+06
 LEFT BOUNDARY OF THE INTERVAL   -9.48134e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.99698e+06  1.60422e+07 |  1.99698e+06  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.99698e+06  1.60422e+07 |  1.99698e+06  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.99698e+06
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           9.93751e+05
  
 NUMBER OF NEGATIVE EIGENVALUES    24
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   9.93751e+05
 LEFT BOUNDARY OF THE INTERVAL   -9.48134e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  9.93751e+05  1.60422e+07 |  9.93751e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  9.93751e+05  1.60422e+07 |  9.93751e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         9.93751e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           4.92135e+05
  
 NUMBER OF NEGATIVE EIGENVALUES    24
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   4.92135e+05
 LEFT BOUNDARY OF THE INTERVAL   -9.48134e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  4.92135e+05  1.60422e+07 |  4.92135e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  4.92135e+05  1.60422e+07 |  4.92135e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.92135e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           2.41327e+05
  
 NUMBER OF NEGATIVE EIGENVALUES    24
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   2.41327e+05
 LEFT BOUNDARY OF THE INTERVAL   -9.48134e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  2.41327e+05  1.60422e+07 |  2.41327e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  2.41327e+05  1.60422e+07 |  2.41327e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         2.41327e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.15923e+05
  
 NUMBER OF NEGATIVE EIGENVALUES    24
 NUMBER OF POSITIVE EIGENVALUES    12
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   1.15923e+05
 LEFT BOUNDARY OF THE INTERVAL   -9.48134e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.15923e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           5.32207e+04
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     8
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   1.15923e+05
 LEFT BOUNDARY OF THE INTERVAL   -9.48134e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  5.32207e+04  5.32207e+04 |  5.32207e+04  5.32207e+04 | -9.48134e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  5.32207e+04  5.32207e+04 |  5.32207e+04  5.32207e+04 | -9.48134e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.41974e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           2.18697e+04
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   2.18697e+04
 LEFT BOUNDARY OF THE INTERVAL   -9.48134e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  2.18697e+04  5.32207e+04 |  2.18697e+04  5.32207e+04 | -9.48134e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  2.18697e+04  5.32207e+04 |  2.18697e+04  5.32207e+04 | -9.48134e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         2.18697e+04
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           6.19416e+03
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   6.19416e+03
 LEFT BOUNDARY OF THE INTERVAL   -9.48134e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  6.19416e+03  5.32207e+04 |  6.19416e+03  5.32207e+04 | -9.48134e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -9.48134e+03 | -3.95574e+37 -9.48134e+03 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  6.19416e+03  5.32207e+04 |  6.19416e+03  5.32207e+04 | -9.48134e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.19416e+03
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -1.64359e+03
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    30
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   6.19416e+03
 LEFT BOUNDARY OF THE INTERVAL   -1.64359e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  6.19416e+03  5.32207e+04 |  6.19416e+03  5.32207e+04 | -9.48134e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.64359e+03 | -3.95574e+37 -1.64359e+03 | -3.95574e+37  6.19416e+03 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  6.19416e+03  5.32207e+04 |  6.19416e+03  5.32207e+04 | -9.48134e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.64359e+03 | -3.95574e+37 -1.64359e+03 | -3.95574e+37  6.19416e+03 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.41974e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -1.64359e+03
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           2.27528e+03
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   2.27528e+03
 LEFT BOUNDARY OF THE INTERVAL   -1.64359e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.64359e+03 | -3.95574e+37 -1.64359e+03 | -3.95574e+37  6.19416e+03 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  2.27528e+03  5.32207e+04 |  2.27528e+03  5.32207e+04 | -1.64359e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.64359e+03 | -3.95574e+37 -1.64359e+03 | -3.95574e+37  6.19416e+03 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  2.27528e+03  5.32207e+04 |  2.27528e+03  5.32207e+04 | -1.64359e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         2.27528e+03
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           3.15845e+02
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   3.15845e+02
 LEFT BOUNDARY OF THE INTERVAL   -1.64359e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.64359e+03 | -3.95574e+37 -1.64359e+03 | -3.95574e+37  6.19416e+03 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  3.15845e+02  5.32207e+04 |  3.15845e+02  5.32207e+04 | -1.64359e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.64359e+03 | -3.95574e+37 -1.64359e+03 | -3.95574e+37  6.19416e+03 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  3.15845e+02  5.32207e+04 |  3.15845e+02  5.32207e+04 | -1.64359e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         3.15845e+02
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -6.63874e+02
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    30
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   3.15845e+02
 LEFT BOUNDARY OF THE INTERVAL   -6.63874e+02
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  3.15845e+02  5.32207e+04 |  3.15845e+02  5.32207e+04 | -1.64359e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -6.63874e+02 | -3.95574e+37 -6.63874e+02 | -3.95574e+37  3.15845e+02 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  3.15845e+02  5.32207e+04 |  3.15845e+02  5.32207e+04 | -1.64359e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -6.63874e+02 | -3.95574e+37 -6.63874e+02 | -3.95574e+37  3.15845e+02 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.41974e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.63874e+02
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -1.74015e+02
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    30
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   3.15845e+02
 LEFT BOUNDARY OF THE INTERVAL   -1.74015e+02
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  3.15845e+02  5.32207e+04 |  3.15845e+02  5.32207e+04 | -1.64359e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.74015e+02 | -3.95574e+37 -1.74015e+02 | -3.95574e+37  3.15845e+02 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  3.15845e+02  5.32207e+04 |  3.15845e+02  5.32207e+04 | -1.64359e+03  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.74015e+02 | -3.95574e+37 -1.74015e+02 | -3.95574e+37  3.15845e+02 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.41974e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -1.74015e+02
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           7.09149e+01
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   7.09149e+01
 LEFT BOUNDARY OF THE INTERVAL   -1.74015e+02
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.74015e+02 | -3.95574e+37 -1.74015e+02 | -3.95574e+37  3.15845e+02 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  7.09149e+01  5.32207e+04 |  7.09149e+01  5.32207e+04 | -1.74015e+02  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.74015e+02 | -3.95574e+37 -1.74015e+02 | -3.95574e+37  3.15845e+02 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  7.09149e+01  5.32207e+04 |  7.09149e+01  5.32207e+04 | -1.74015e+02  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         7.09149e+01
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -5.15500e+01
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    30
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   7.09149e+01
 LEFT BOUNDARY OF THE INTERVAL   -5.15500e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  7.09149e+01  5.32207e+04 |  7.09149e+01  5.32207e+04 | -1.74015e+02  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -5.15500e+01 | -3.95574e+37 -5.15500e+01 | -3.95574e+37  7.09149e+01 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  7.09149e+01  5.32207e+04 |  7.09149e+01  5.32207e+04 | -1.74015e+02  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -5.15500e+01 | -3.95574e+37 -5.15500e+01 | -3.95574e+37  7.09149e+01 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.41974e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -5.15500e+01
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           9.68248e+00
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   9.68248e+00
 LEFT BOUNDARY OF THE INTERVAL   -5.15500e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -5.15500e+01 | -3.95574e+37 -5.15500e+01 | -3.95574e+37  7.09149e+01 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  9.68248e+00  5.32207e+04 |  9.68248e+00  5.32207e+04 | -5.15500e+01  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -5.15500e+01 | -3.95574e+37 -5.15500e+01 | -3.95574e+37  7.09149e+01 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  9.68248e+00  5.32207e+04 |  9.68248e+00  5.32207e+04 | -5.15500e+01  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         9.68248e+00
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -2.09337e+01
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    30
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   9.68248e+00
 LEFT BOUNDARY OF THE INTERVAL   -2.09337e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  9.68248e+00  5.32207e+04 |  9.68248e+00  5.32207e+04 | -5.15500e+01  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.09337e+01 | -3.95574e+37 -2.09337e+01 | -3.95574e+37  9.68248e+00 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  9.68248e+00  5.32207e+04 |  9.68248e+00  5.32207e+04 | -5.15500e+01  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.09337e+01 | -3.95574e+37 -2.09337e+01 | -3.95574e+37  9.68248e+00 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.41974e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -2.09337e+01
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -5.62563e+00
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    30
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   9.68248e+00
 LEFT BOUNDARY OF THE INTERVAL   -5.62563e+00
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  9.68248e+00  5.32207e+04 |  9.68248e+00  5.32207e+04 | -5.15500e+01  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -5.62563e+00 | -3.95574e+37 -5.62563e+00 | -3.95574e+37  9.68248e+00 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  9.68248e+00  5.32207e+04 |  9.68248e+00  5.32207e+04 | -5.15500e+01  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -5.62563e+00 | -3.95574e+37 -5.62563e+00 | -3.95574e+37  9.68248e+00 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.41974e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -5.62563e+00
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           2.02843e+00
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   2.02843e+00
 LEFT BOUNDARY OF THE INTERVAL   -5.62563e+00
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -5.62563e+00 | -3.95574e+37 -5.62563e+00 | -3.95574e+37  9.68248e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  2.02843e+00  5.32207e+04 |  2.02843e+00  5.32207e+04 | -5.62563e+00  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -5.62563e+00 | -3.95574e+37 -5.62563e+00 | -3.95574e+37  9.68248e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  2.02843e+00  5.32207e+04 |  2.02843e+00  5.32207e+04 | -5.62563e+00  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         2.02843e+00
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -1.79860e+00
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    30
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   2.02843e+00
 LEFT BOUNDARY OF THE INTERVAL   -1.79860e+00
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  2.02843e+00  5.32207e+04 |  2.02843e+00  5.32207e+04 | -5.62563e+00  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.79860e+00 | -3.95574e+37 -1.79860e+00 | -3.95574e+37  2.02843e+00 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  2.02843e+00  5.32207e+04 |  2.02843e+00  5.32207e+04 | -5.62563e+00  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.79860e+00 | -3.95574e+37 -1.79860e+00 | -3.95574e+37  2.02843e+00 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.41974e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -1.79860e+00
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.14913e-01
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   1.14913e-01
 LEFT BOUNDARY OF THE INTERVAL   -1.79860e+00
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.79860e+00 | -3.95574e+37 -1.79860e+00 | -3.95574e+37  2.02843e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  1.14913e-01  5.32207e+04 |  1.14913e-01  5.32207e+04 | -1.79860e+00  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.79860e+00 | -3.95574e+37 -1.79860e+00 | -3.95574e+37  2.02843e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  1.14913e-01  5.32207e+04 |  1.14913e-01  5.32207e+04 | -1.79860e+00  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.14913e-01
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -8.41844e-01
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    30
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.14913e-01
 LEFT BOUNDARY OF THE INTERVAL   -8.41844e-01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  1.14913e-01  5.32207e+04 |  1.14913e-01  5.32207e+04 | -1.79860e+00  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -8.41844e-01 | -3.95574e+37 -8.41844e-01 | -3.95574e+37  1.14913e-01 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  1.14913e-01  5.32207e+04 |  1.14913e-01  5.32207e+04 | -1.79860e+00  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -8.41844e-01 | -3.95574e+37 -8.41844e-01 | -3.95574e+37  1.14913e-01 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.41974e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -8.41844e-01
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -3.63466e-01
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    30
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.14913e-01
 LEFT BOUNDARY OF THE INTERVAL   -3.63466e-01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  1.14913e-01  5.32207e+04 |  1.14913e-01  5.32207e+04 | -1.79860e+00  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -3.63466e-01 | -3.95574e+37 -3.63466e-01 | -3.95574e+37  1.14913e-01 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  1.14913e-01  5.32207e+04 |  1.14913e-01  5.32207e+04 | -1.79860e+00  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -3.63466e-01 | -3.95574e+37 -3.63466e-01 | -3.95574e+37  1.14913e-01 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.41974e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -3.63466e-01
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -1.24276e-01
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    30
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.14913e-01
 LEFT BOUNDARY OF THE INTERVAL   -1.24276e-01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  1.14913e-01  5.32207e+04 |  1.14913e-01  5.32207e+04 | -1.79860e+00  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.24276e-01 | -3.95574e+37 -1.24276e-01 | -3.95574e+37  1.14913e-01 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  1.14913e-01  5.32207e+04 |  1.14913e-01  5.32207e+04 | -1.79860e+00  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.24276e-01 | -3.95574e+37 -1.24276e-01 | -3.95574e+37  1.14913e-01 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.41974e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -1.24276e-01
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -4.68180e-03
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    30
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.14913e-01
 LEFT BOUNDARY OF THE INTERVAL   -4.68180e-03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       12       12 |     0     0 |     0     0 |      6     8       |
 |  1.14913e-01  5.32207e+04 |  1.14913e-01  5.32207e+04 | -1.79860e+00  1.15923e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -4.68180e-03 | -3.95574e+37 -4.68180e-03 | -3.95574e+37  1.14913e-01 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   10
 ALL THE EIGENVALUES REQUESTED ON THIS LANCZOS RUN ARE FOUND
  
 NUMBER OF ACCEPTED EIGENVALUES   14
  
 ACCEPTED EIGENVALUES
  -1.46085e-09 -9.60703e-10 -7.48960e-10 -8.64476e-11  3.14462e-10  1.46730e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     0 |      0     8       |
 | -3.95574e+37  5.32207e+04 | -3.95574e+37  5.32207e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
 SENTINEL WAS NOT INCLUDED INTO THE TRUST INTERVAL LIST. SHIFT  -4.68180e-03
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    6
 SHIFT RIGHT         1.88679e+04
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.68180e-03
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           7.85461e+04
  
 NUMBER OF NEGATIVE EIGENVALUES    12
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     8
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.15923e+05
 LEFT BOUNDARY OF THE INTERVAL    7.85461e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     0 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   15
  
 ACCEPTED EIGENVALUES
  -1.46085e-09 -9.60703e-10 -7.48960e-10 -8.64476e-11  3.14462e-10  1.46730e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07
   1.97239e+07  1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     0 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.07266e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.85461e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           7.85461e+04
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     7
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.15923e+05
 LEFT BOUNDARY OF THE INTERVAL    7.85461e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   16
  
 ACCEPTED EIGENVALUES
  -1.46085e-09 -9.60703e-10 -7.48960e-10 -8.64476e-11  3.14462e-10  1.46730e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   1.97239e+07  1.97239e+07  1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.07266e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.85461e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           7.85461e+04
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.15923e+05
 LEFT BOUNDARY OF THE INTERVAL    7.85461e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   17
  
 ACCEPTED EIGENVALUES
  -1.46085e-09 -9.60703e-10 -7.48960e-10 -8.64476e-11  3.14462e-10  1.46730e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.29783e+04  1.97239e+07  1.97239e+07  1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.07266e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.85461e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           7.85461e+04
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     5
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.15923e+05
 LEFT BOUNDARY OF THE INTERVAL    7.85461e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   18
  
 ACCEPTED EIGENVALUES
  -1.46085e-09 -9.60703e-10 -7.48960e-10 -8.64476e-11  3.14462e-10  1.46730e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07  1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.07266e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.85461e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           7.85461e+04
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     4
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.15923e+05
 LEFT BOUNDARY OF THE INTERVAL    7.85461e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   19
  
 ACCEPTED EIGENVALUES
  -1.46085e-09 -9.60703e-10 -7.48960e-10 -8.64476e-11  3.14462e-10  1.46730e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07  1.97239e+07
   6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.07266e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.85461e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           7.85461e+04
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     3
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.15923e+05
 LEFT BOUNDARY OF THE INTERVAL    7.85461e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   20
  
 ACCEPTED EIGENVALUES
  -1.46085e-09 -9.60703e-10 -7.48960e-10 -8.64476e-11  3.14462e-10  1.46730e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.07266e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.85461e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           7.85461e+04
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     2
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.15923e+05
 LEFT BOUNDARY OF THE INTERVAL    7.85461e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   21
  
 ACCEPTED EIGENVALUES
  -1.46085e-09 -9.60703e-10 -7.48960e-10 -8.64476e-11  3.14462e-10  1.46730e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07
   1.97239e+07  1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.07266e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.85461e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           7.85461e+04
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     1
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.15923e+05
 LEFT BOUNDARY OF THE INTERVAL    7.85461e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     14     2       |
 |  1.15923e+05  1.60422e+07 |  1.15923e+05  1.60422e+07 | -9.48134e+03  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        6       12 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.85461e+04 | -3.95574e+37  7.85461e+04 | -3.95574e+37  1.15923e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES REQUESTED ON THIS LANCZOS RUN ARE FOUND
  
 NUMBER OF ACCEPTED EIGENVALUES   22
  
 ACCEPTED EIGENVALUES
  -1.46085e-09 -9.60703e-10 -7.48960e-10 -8.64476e-11  3.14462e-10  1.46730e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   1.97239e+07  1.97239e+07  1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    2
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |     16     0       |
 |  3.20939e+07  6.41974e+07 |  3.20939e+07  6.41974e+07 | -9.48134e+03  6.41974e+07 |
 -------------------------------------------------------------------------------------
 |  +  |     18 |        6       24 |     0     0 |     0     0 |      0     2       |
 | -3.95574e+37  1.60422e+07 | -3.95574e+37  1.60422e+07 | -3.95574e+37  3.20939e+07 |
 -------------------------------------------------------------------------------------
 SENTINEL WAS NOT INCLUDED INTO THE TRUST INTERVAL LIST. SHIFT   7.85461e+04
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL    6
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.07266e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.85461e+04
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           2.26905e+07
  
 NUMBER OF NEGATIVE EIGENVALUES    29
 NUMBER OF POSITIVE EIGENVALUES     7
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    23
 NUMBER OF LAGRANGE MULTIPLIERS     6
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      2
 RIGHT BOUNDARY OF THE INTERVAL   2.26905e+07
 LEFT BOUNDARY OF THE INTERVAL    1.60422e+07
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |     18 |        6       24 |     0     0 |     0     0 |      0     2       |
 | -3.95574e+37  1.60422e+07 | -3.95574e+37  1.60422e+07 | -3.95574e+37  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |      2     0       |
 |  2.26905e+07  6.41974e+07 |  2.26905e+07  6.41974e+07 |  1.60422e+07  6.41974e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    4
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   23
  
 ACCEPTED EIGENVALUES
  -1.46085e-09 -9.60703e-10 -7.48960e-10 -8.64476e-11  3.14462e-10  1.46730e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   1.97239e+07  1.97239e+07  1.97239e+07  1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    2
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |     18 |        6       24 |     0     0 |     0     0 |      0     2       |
 | -3.95574e+37  1.60422e+07 | -3.95574e+37  1.60422e+07 | -3.95574e+37  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     0     0 |      2     0       |
 |  2.26905e+07  6.41974e+07 |  2.26905e+07  6.41974e+07 |  1.60422e+07  6.41974e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL    6
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         2.26905e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     1
 SHIFT LEFT          1.67869e+07
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           2.26905e+07
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      1
 RIGHT BOUNDARY OF THE INTERVAL   2.26905e+07
 LEFT BOUNDARY OF THE INTERVAL    1.60422e+07
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |     18 |        6       24 |     0     0 |     0     0 |      0     2       |
 | -3.95574e+37  1.60422e+07 | -3.95574e+37  1.60422e+07 | -3.95574e+37  3.20939e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       29       30 |     0     0 |     1     0 |      2     0       |
 |  2.26905e+07  6.41974e+07 |  2.26905e+07  6.41974e+07 |  1.60422e+07  6.41974e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    4
 ALL THE EIGENVALUES REQUESTED ON THIS LANCZOS RUN ARE FOUND
  
 NUMBER OF ACCEPTED EIGENVALUES   24
  
 ACCEPTED EIGENVALUES
  -1.46085e-09 -9.60703e-10 -7.48960e-10 -8.64476e-11  3.14462e-10  1.46730e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   1.97239e+07  1.97239e+07  1.97239e+07  1.97239e+07  1.97239e+07  6.41026e+07
  
 THE REQUIRED NUMBER OF EIGENVALUES IS FOUND
  
 FINAL TRUST INTERVALS AFTER SUCCESSFUL LANCZOS RUN
  
 NUMBER OF TRUST INTERVALS    1
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |     24 |        6       30 |     0     0 |     0     0 |      0     0       |
 | -3.95574e+37  6.41974e+07 | -3.95574e+37  6.41974e+07 | -3.95574e+37  6.41974e+07 |
 -------------------------------------------------------------------------------------
  
 FIRST EIGENVALUE REQUESTED BY THE USER     1
  
 EIGENVALUES REQUESTED BY THE USER
  -1.46085e-09 -9.60703e-10 -7.48960e-10 -8.64476e-11  3.14462e-10  1.46730e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   1.97239e+07  1.97239e+07  1.97239e+07  1.97239e+07  1.97239e+07  6.41026e+07
  
 EIGENVALUE EXTRACTION USING THE LANCZOS EIGENSOLVER IS COMPLETE
 NUMBER OF LANCZOS RUNS       :  41
 NUMBER OF LANCZOS STEPS      :  84
 NUMBER OF FACTORIZATIONS     :  33
 NUMBER OF BACK SUBSTITUTIONS :  95

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.41000    
       ELAPSED SYSTEM TIME (SEC)    =  8.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.49000    
       ELAPSED WALLCLOCK TIME (SEC) =          1


          THE ANALYSIS HAS BEEN COMPLETED



     ANALYSIS SUMMARY:
     TOTAL OF          2  INCREMENTS
                       0  CUTBACKS IN AUTOMATIC INCREMENTATION
                       0  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                     128  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                      33  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       1  WARNING MESSAGES DURING USER INPUT PROCESSING
                       1  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                       0  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       0  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =   1.2600    
       SYSTEM TIME (SEC)    =  0.18000    
       TOTAL CPU TIME (SEC) =   1.4400    
       WALLCLOCK TIME (SEC) =          4
