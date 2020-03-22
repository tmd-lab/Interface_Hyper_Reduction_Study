1

   Abaqus 3DEXPERIENCE R2019x                                  Date 18-Mar-2020   Time 10:21:00
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

     CHARACTERISTIC ELEMENT LENGTH      2.00    

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     CHECK POINT   START OF SOLVER    
  
 COMPUTER PRECISION USED BY ABAQUS       2.22045e-16
 CONVERGENCE CRITERION FOR LANCZOS RUN   1.81899e-12
  
 NUMBER OF EQUATIONS                                                          24

 ***WARNING: ONLY 24 EIGENMODES WILL BE EXTRACTED BECAUSE THERE ARE NO MORE 
             UNRESTRAINED DEGREES OF FREEDOM
  
 NUMBER OF UNRESTRAINED DEGREES OF FREEDOM                                    24
 HEURISTIC SCALE FOR THE FIRST NONZERO EIGENVALUES (PROBLEM SCALE)   1.18517e+04
 INITIAL SHIFT VALUE                                                -1.18517e+04
  
 NUMBER OF CPUS USED BY DIRECT SOLVER     1
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
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
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.42211e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -1.18517e+04
  
 SHIFTING TO THE RIGHT DIRECTION
 NEW SHIFT           6.42211e+07
  
 NUMBER OF NEGATIVE EIGENVALUES    24
 NUMBER OF POSITIVE EIGENVALUES     0
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    24
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: EXPAND THE TRUST INTERVAL TO THE RIGHT
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     16
 RIGHT BOUNDARY OF THE INTERVAL   6.42211e+07
 LEFT BOUNDARY OF THE INTERVAL   -1.18517e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     16     0       |
 |  6.42211e+07  6.42211e+07 |  6.42211e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
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
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       24       24 |     0     0 |     0     0 |     16     0       |
 |  6.42211e+07  6.42211e+07 |  6.42211e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.42211e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           3.21046e+07
  
 NUMBER OF NEGATIVE EIGENVALUES    23
 NUMBER OF POSITIVE EIGENVALUES     1
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    23
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     16
 RIGHT BOUNDARY OF THE INTERVAL   3.21046e+07
 LEFT BOUNDARY OF THE INTERVAL   -1.18517e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
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
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         3.21046e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.60464e+07
  
 NUMBER OF NEGATIVE EIGENVALUES    18
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     2
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   3.21046e+07
 LEFT BOUNDARY OF THE INTERVAL   -1.18517e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.60464e+07  1.60464e+07 |  1.60464e+07  1.60464e+07 | -1.18517e+04  3.21046e+07 |
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
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.60464e+07  1.60464e+07 |  1.60464e+07  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.42211e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           8.01726e+06
  
 NUMBER OF NEGATIVE EIGENVALUES    18
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   8.01726e+06
 LEFT BOUNDARY OF THE INTERVAL   -1.18517e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  8.01726e+06  1.60464e+07 |  8.01726e+06  1.60464e+07 | -1.18517e+04  3.21046e+07 |
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
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  8.01726e+06  1.60464e+07 |  8.01726e+06  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         8.01726e+06
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           4.00271e+06
  
 NUMBER OF NEGATIVE EIGENVALUES    18
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   4.00271e+06
 LEFT BOUNDARY OF THE INTERVAL   -1.18517e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  4.00271e+06  1.60464e+07 |  4.00271e+06  1.60464e+07 | -1.18517e+04  3.21046e+07 |
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
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  4.00271e+06  1.60464e+07 |  4.00271e+06  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.00271e+06
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.99543e+06
  
 NUMBER OF NEGATIVE EIGENVALUES    18
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   1.99543e+06
 LEFT BOUNDARY OF THE INTERVAL   -1.18517e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.99543e+06  1.60464e+07 |  1.99543e+06  1.60464e+07 | -1.18517e+04  3.21046e+07 |
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
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.99543e+06  1.60464e+07 |  1.99543e+06  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.99543e+06
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           9.91788e+05
  
 NUMBER OF NEGATIVE EIGENVALUES    18
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   9.91788e+05
 LEFT BOUNDARY OF THE INTERVAL   -1.18517e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  9.91788e+05  1.60464e+07 |  9.91788e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
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
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  9.91788e+05  1.60464e+07 |  9.91788e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         9.91788e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           4.89968e+05
  
 NUMBER OF NEGATIVE EIGENVALUES    18
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   4.89968e+05
 LEFT BOUNDARY OF THE INTERVAL   -1.18517e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  4.89968e+05  1.60464e+07 |  4.89968e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
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
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  4.89968e+05  1.60464e+07 |  4.89968e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         4.89968e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           2.39058e+05
  
 NUMBER OF NEGATIVE EIGENVALUES    18
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   2.39058e+05
 LEFT BOUNDARY OF THE INTERVAL   -1.18517e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  2.39058e+05  1.60464e+07 |  2.39058e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
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
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  2.39058e+05  1.60464e+07 |  2.39058e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         2.39058e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.13603e+05
  
 NUMBER OF NEGATIVE EIGENVALUES    18
 NUMBER OF POSITIVE EIGENVALUES     6
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    18
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     14
 RIGHT BOUNDARY OF THE INTERVAL   1.13603e+05
 LEFT BOUNDARY OF THE INTERVAL   -1.18517e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
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
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.13603e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           5.08758e+04
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     8
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   1.13603e+05
 LEFT BOUNDARY OF THE INTERVAL   -1.18517e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  5.08758e+04  5.08758e+04 |  5.08758e+04  5.08758e+04 | -1.18517e+04  1.13603e+05 |
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
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  5.08758e+04  5.08758e+04 |  5.08758e+04  5.08758e+04 | -1.18517e+04  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.42211e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.95121e+04
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   1.95121e+04
 LEFT BOUNDARY OF THE INTERVAL   -1.18517e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.95121e+04  5.08758e+04 |  1.95121e+04  5.08758e+04 | -1.18517e+04  1.13603e+05 |
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
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.95121e+04  5.08758e+04 |  1.95121e+04  5.08758e+04 | -1.18517e+04  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.95121e+04
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           3.83019e+03
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   3.83019e+03
 LEFT BOUNDARY OF THE INTERVAL   -1.18517e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  3.83019e+03  5.08758e+04 |  3.83019e+03  5.08758e+04 | -1.18517e+04  1.13603e+05 |
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
 |  -  |      0 |        0        0 |     0     0 |     0     8 |      0    24       |
 | -3.95574e+37 -1.18517e+04 | -3.95574e+37 -1.18517e+04 | -3.95574e+37  3.95574e+37 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  3.83019e+03  5.08758e+04 |  3.83019e+03  5.08758e+04 | -1.18517e+04  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         3.83019e+03
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -4.01075e+03
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   3.83019e+03
 LEFT BOUNDARY OF THE INTERVAL   -4.01075e+03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  3.83019e+03  5.08758e+04 |  3.83019e+03  5.08758e+04 | -1.18517e+04  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -4.01075e+03 | -3.95574e+37 -4.01075e+03 | -3.95574e+37  3.83019e+03 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  3.83019e+03  5.08758e+04 |  3.83019e+03  5.08758e+04 | -1.18517e+04  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -4.01075e+03 | -3.95574e+37 -4.01075e+03 | -3.95574e+37  3.83019e+03 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.42211e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -4.01075e+03
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -9.02794e+01
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   3.83019e+03
 LEFT BOUNDARY OF THE INTERVAL   -9.02794e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  3.83019e+03  5.08758e+04 |  3.83019e+03  5.08758e+04 | -1.18517e+04  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -9.02794e+01 | -3.95574e+37 -9.02794e+01 | -3.95574e+37  3.83019e+03 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  3.83019e+03  5.08758e+04 |  3.83019e+03  5.08758e+04 | -1.18517e+04  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -9.02794e+01 | -3.95574e+37 -9.02794e+01 | -3.95574e+37  3.83019e+03 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.42211e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -9.02794e+01
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.86995e+03
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   1.86995e+03
 LEFT BOUNDARY OF THE INTERVAL   -9.02794e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -9.02794e+01 | -3.95574e+37 -9.02794e+01 | -3.95574e+37  3.83019e+03 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.86995e+03  5.08758e+04 |  1.86995e+03  5.08758e+04 | -9.02794e+01  1.13603e+05 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -9.02794e+01 | -3.95574e+37 -9.02794e+01 | -3.95574e+37  3.83019e+03 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.86995e+03  5.08758e+04 |  1.86995e+03  5.08758e+04 | -9.02794e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.86995e+03
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           8.89837e+02
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   8.89837e+02
 LEFT BOUNDARY OF THE INTERVAL   -9.02794e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -9.02794e+01 | -3.95574e+37 -9.02794e+01 | -3.95574e+37  3.83019e+03 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  8.89837e+02  5.08758e+04 |  8.89837e+02  5.08758e+04 | -9.02794e+01  1.13603e+05 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -9.02794e+01 | -3.95574e+37 -9.02794e+01 | -3.95574e+37  3.83019e+03 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  8.89837e+02  5.08758e+04 |  8.89837e+02  5.08758e+04 | -9.02794e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         8.89837e+02
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           3.99779e+02
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   3.99779e+02
 LEFT BOUNDARY OF THE INTERVAL   -9.02794e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -9.02794e+01 | -3.95574e+37 -9.02794e+01 | -3.95574e+37  3.83019e+03 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  3.99779e+02  5.08758e+04 |  3.99779e+02  5.08758e+04 | -9.02794e+01  1.13603e+05 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -9.02794e+01 | -3.95574e+37 -9.02794e+01 | -3.95574e+37  3.83019e+03 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  3.99779e+02  5.08758e+04 |  3.99779e+02  5.08758e+04 | -9.02794e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         3.99779e+02
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.54750e+02
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   1.54750e+02
 LEFT BOUNDARY OF THE INTERVAL   -9.02794e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -9.02794e+01 | -3.95574e+37 -9.02794e+01 | -3.95574e+37  3.83019e+03 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.54750e+02  5.08758e+04 |  1.54750e+02  5.08758e+04 | -9.02794e+01  1.13603e+05 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -9.02794e+01 | -3.95574e+37 -9.02794e+01 | -3.95574e+37  3.83019e+03 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.54750e+02  5.08758e+04 |  1.54750e+02  5.08758e+04 | -9.02794e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.54750e+02
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           3.22352e+01
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   3.22352e+01
 LEFT BOUNDARY OF THE INTERVAL   -9.02794e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -9.02794e+01 | -3.95574e+37 -9.02794e+01 | -3.95574e+37  3.83019e+03 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  3.22352e+01  5.08758e+04 |  3.22352e+01  5.08758e+04 | -9.02794e+01  1.13603e+05 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -9.02794e+01 | -3.95574e+37 -9.02794e+01 | -3.95574e+37  3.83019e+03 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  3.22352e+01  5.08758e+04 |  3.22352e+01  5.08758e+04 | -9.02794e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         3.22352e+01
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -2.90221e+01
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   3.22352e+01
 LEFT BOUNDARY OF THE INTERVAL   -2.90221e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  3.22352e+01  5.08758e+04 |  3.22352e+01  5.08758e+04 | -9.02794e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.90221e+01 | -3.95574e+37 -2.90221e+01 | -3.95574e+37  3.22352e+01 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  3.22352e+01  5.08758e+04 |  3.22352e+01  5.08758e+04 | -9.02794e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.90221e+01 | -3.95574e+37 -2.90221e+01 | -3.95574e+37  3.22352e+01 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.42211e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -2.90221e+01
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.60654e+00
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   1.60654e+00
 LEFT BOUNDARY OF THE INTERVAL   -2.90221e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.90221e+01 | -3.95574e+37 -2.90221e+01 | -3.95574e+37  3.22352e+01 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.60654e+00  5.08758e+04 |  1.60654e+00  5.08758e+04 | -2.90221e+01  1.13603e+05 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.90221e+01 | -3.95574e+37 -2.90221e+01 | -3.95574e+37  3.22352e+01 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.60654e+00  5.08758e+04 |  1.60654e+00  5.08758e+04 | -2.90221e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.60654e+00
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -1.37078e+01
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.60654e+00
 LEFT BOUNDARY OF THE INTERVAL   -1.37078e+01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.60654e+00  5.08758e+04 |  1.60654e+00  5.08758e+04 | -2.90221e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.37078e+01 | -3.95574e+37 -1.37078e+01 | -3.95574e+37  1.60654e+00 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.60654e+00  5.08758e+04 |  1.60654e+00  5.08758e+04 | -2.90221e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -1.37078e+01 | -3.95574e+37 -1.37078e+01 | -3.95574e+37  1.60654e+00 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.42211e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -1.37078e+01
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -6.05062e+00
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.60654e+00
 LEFT BOUNDARY OF THE INTERVAL   -6.05062e+00
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.60654e+00  5.08758e+04 |  1.60654e+00  5.08758e+04 | -2.90221e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -6.05062e+00 | -3.95574e+37 -6.05062e+00 | -3.95574e+37  1.60654e+00 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.60654e+00  5.08758e+04 |  1.60654e+00  5.08758e+04 | -2.90221e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -6.05062e+00 | -3.95574e+37 -6.05062e+00 | -3.95574e+37  1.60654e+00 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.42211e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.05062e+00
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -2.22204e+00
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.60654e+00
 LEFT BOUNDARY OF THE INTERVAL   -2.22204e+00
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.60654e+00  5.08758e+04 |  1.60654e+00  5.08758e+04 | -2.90221e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.22204e+00 | -3.95574e+37 -2.22204e+00 | -3.95574e+37  1.60654e+00 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.60654e+00  5.08758e+04 |  1.60654e+00  5.08758e+04 | -2.90221e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -2.22204e+00 | -3.95574e+37 -2.22204e+00 | -3.95574e+37  1.60654e+00 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.42211e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -2.22204e+00
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -3.07750e-01
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.60654e+00
 LEFT BOUNDARY OF THE INTERVAL   -3.07750e-01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.60654e+00  5.08758e+04 |  1.60654e+00  5.08758e+04 | -2.90221e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -3.07750e-01 | -3.95574e+37 -3.07750e-01 | -3.95574e+37  1.60654e+00 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.60654e+00  5.08758e+04 |  1.60654e+00  5.08758e+04 | -2.90221e+01  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -3.07750e-01 | -3.95574e+37 -3.07750e-01 | -3.95574e+37  1.60654e+00 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.42211e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -3.07750e-01
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           6.49395e-01
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   6.49395e-01
 LEFT BOUNDARY OF THE INTERVAL   -3.07750e-01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -3.07750e-01 | -3.95574e+37 -3.07750e-01 | -3.95574e+37  1.60654e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  6.49395e-01  5.08758e+04 |  6.49395e-01  5.08758e+04 | -3.07750e-01  1.13603e+05 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -3.07750e-01 | -3.95574e+37 -3.07750e-01 | -3.95574e+37  1.60654e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  6.49395e-01  5.08758e+04 |  6.49395e-01  5.08758e+04 | -3.07750e-01  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.49395e-01
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           1.70822e-01
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   1.70822e-01
 LEFT BOUNDARY OF THE INTERVAL   -3.07750e-01
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -3.07750e-01 | -3.95574e+37 -3.07750e-01 | -3.95574e+37  1.60654e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.70822e-01  5.08758e+04 |  1.70822e-01  5.08758e+04 | -3.07750e-01  1.13603e+05 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -3.07750e-01 | -3.95574e+37 -3.07750e-01 | -3.95574e+37  1.60654e+00 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.70822e-01  5.08758e+04 |  1.70822e-01  5.08758e+04 | -3.07750e-01  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         1.70822e-01
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -6.84639e-02
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.70822e-01
 LEFT BOUNDARY OF THE INTERVAL   -6.84639e-02
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.70822e-01  5.08758e+04 |  1.70822e-01  5.08758e+04 | -3.07750e-01  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -6.84639e-02 | -3.95574e+37 -6.84639e-02 | -3.95574e+37  1.70822e-01 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  1.70822e-01  5.08758e+04 |  1.70822e-01  5.08758e+04 | -3.07750e-01  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -6.84639e-02 | -3.95574e+37 -6.84639e-02 | -3.95574e+37  1.70822e-01 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         6.42211e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.84639e-02
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           5.11793e-02
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      6
 RIGHT BOUNDARY OF THE INTERVAL   5.11793e-02
 LEFT BOUNDARY OF THE INTERVAL   -6.84639e-02
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -6.84639e-02 | -3.95574e+37 -6.84639e-02 | -3.95574e+37  1.70822e-01 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  5.11793e-02  5.08758e+04 |  5.11793e-02  5.08758e+04 | -6.84639e-02  1.13603e+05 |
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -6.84639e-02 | -3.95574e+37 -6.84639e-02 | -3.95574e+37  1.70822e-01 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  5.11793e-02  5.08758e+04 |  5.11793e-02  5.08758e+04 | -6.84639e-02  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   24
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         5.11793e-02
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -6.41026e+07
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          -8.64227e-03
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES    24
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   5.11793e-02
 LEFT BOUNDARY OF THE INTERVAL   -8.64227e-03
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        6        6 |     0     0 |     0     0 |      6     8       |
 |  5.11793e-02  5.08758e+04 |  5.11793e-02  5.08758e+04 | -6.84639e-02  1.13603e+05 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |        0        0 |     0     0 |     0     0 |      0     6       |
 | -3.95574e+37 -8.64227e-03 | -3.95574e+37 -8.64227e-03 | -3.95574e+37  5.11793e-02 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   10
 ALL THE EIGENVALUES REQUESTED ON THIS LANCZOS RUN ARE FOUND
  
 NUMBER OF ACCEPTED EIGENVALUES   14
  
 ACCEPTED EIGENVALUES
  -1.16568e-09 -6.58454e-10 -3.27990e-10 -6.41666e-11  8.25525e-10  1.37336e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07
   1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     0 |      0     8       |
 | -3.95574e+37  5.08758e+04 | -3.95574e+37  5.08758e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
 SENTINEL WAS NOT INCLUDED INTO THE TRUST INTERVAL LIST. SHIFT  -8.64227e-03
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    6
 SHIFT RIGHT         2.35848e+04
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         -8.64227e-03
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           7.60241e+04
  
 NUMBER OF NEGATIVE EIGENVALUES     6
 NUMBER OF POSITIVE EIGENVALUES    18
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     6
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     8
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.13603e+05
 LEFT BOUNDARY OF THE INTERVAL    7.60241e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     0 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   15
  
 ACCEPTED EIGENVALUES
  -1.16568e-09 -6.58454e-10 -3.27990e-10 -6.41666e-11  8.25525e-10  1.37336e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  1.97239e+07
   1.97239e+07  1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     0 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.09763e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.60241e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           7.60241e+04
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     7
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.13603e+05
 LEFT BOUNDARY OF THE INTERVAL    7.60241e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   16
  
 ACCEPTED EIGENVALUES
  -1.16568e-09 -6.58454e-10 -3.27990e-10 -6.41666e-11  8.25525e-10  1.37336e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   1.97239e+07  1.97239e+07  1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.09763e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.60241e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           7.60241e+04
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     6
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.13603e+05
 LEFT BOUNDARY OF THE INTERVAL    7.60241e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   17
  
 ACCEPTED EIGENVALUES
  -1.16568e-09 -6.58454e-10 -3.27990e-10 -6.41666e-11  8.25525e-10  1.37336e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.29783e+04  1.97239e+07  1.97239e+07  1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.09763e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.60241e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           7.60241e+04
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     5
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.13603e+05
 LEFT BOUNDARY OF THE INTERVAL    7.60241e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   18
  
 ACCEPTED EIGENVALUES
  -1.16568e-09 -6.58454e-10 -3.27990e-10 -6.41666e-11  8.25525e-10  1.37336e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.29783e+04  9.29783e+04  1.97239e+07  1.97239e+07  1.97239e+07  6.41026e+07
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.09763e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.60241e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           7.60241e+04
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     4
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.13603e+05
 LEFT BOUNDARY OF THE INTERVAL    7.60241e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   19
  
 ACCEPTED EIGENVALUES
  -1.16568e-09 -6.58454e-10 -3.27990e-10 -6.41666e-11  8.25525e-10  1.37336e-09
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.09763e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.60241e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           7.60241e+04
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     3
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.13603e+05
 LEFT BOUNDARY OF THE INTERVAL    7.60241e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   20
  
 ACCEPTED EIGENVALUES
  -1.16568e-09 -6.58454e-10 -3.27990e-10 -6.41666e-11  8.25525e-10  1.37336e-09
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.09763e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.60241e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           7.60241e+04
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     2
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.13603e+05
 LEFT BOUNDARY OF THE INTERVAL    7.60241e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   21
  
 ACCEPTED EIGENVALUES
  -1.16568e-09 -6.58454e-10 -3.27990e-10 -6.41666e-11  8.25525e-10  1.37336e-09
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   18
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.09763e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.60241e+04
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           7.60241e+04
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     1
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL   1.13603e+05
 LEFT BOUNDARY OF THE INTERVAL    7.60241e+04
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       18       18 |     0     0 |     0     0 |     14     2       |
 |  1.13603e+05  1.60464e+07 |  1.13603e+05  1.60464e+07 | -1.18517e+04  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  +  |      6 |        0        6 |     0     0 |     0     1 |      0     8       |
 | -3.95574e+37  7.60241e+04 | -3.95574e+37  7.60241e+04 | -3.95574e+37  1.13603e+05 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    3
 ALL THE EIGENVALUES REQUESTED ON THIS LANCZOS RUN ARE FOUND
  
 NUMBER OF ACCEPTED EIGENVALUES   22
  
 ACCEPTED EIGENVALUES
  -1.16568e-09 -6.58454e-10 -3.27990e-10 -6.41666e-11  8.25525e-10  1.37336e-09
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
 |  -  |      1 |       23       24 |     0     0 |     0     0 |     16     0       |
 |  3.21046e+07  6.42211e+07 |  3.21046e+07  6.42211e+07 | -1.18517e+04  6.42211e+07 |
 -------------------------------------------------------------------------------------
 |  +  |     18 |        0       18 |     0     0 |     0     0 |      0     2       |
 | -3.95574e+37  1.60464e+07 | -3.95574e+37  1.60464e+07 | -3.95574e+37  3.21046e+07 |
 -------------------------------------------------------------------------------------
 SENTINEL WAS NOT INCLUDED INTO THE TRUST INTERVAL LIST. SHIFT   7.60241e+04
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL    6
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    1
 SHIFT RIGHT         1.09763e+05
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT          7.60241e+04
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT           2.26972e+07
  
 NUMBER OF NEGATIVE EIGENVALUES    23
 NUMBER OF POSITIVE EIGENVALUES     1
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    23
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      2
 RIGHT BOUNDARY OF THE INTERVAL   2.26972e+07
 LEFT BOUNDARY OF THE INTERVAL    1.60464e+07
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |     18 |        0       18 |     0     0 |     0     0 |      0     2       |
 | -3.95574e+37  1.60464e+07 | -3.95574e+37  1.60464e+07 | -3.95574e+37  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |      2     0       |
 |  2.26972e+07  6.42211e+07 |  2.26972e+07  6.42211e+07 |  1.60464e+07  6.42211e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    4
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   23
  
 ACCEPTED EIGENVALUES
  -1.16568e-09 -6.58454e-10 -3.27990e-10 -6.41666e-11  8.25525e-10  1.37336e-09
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
 |  +  |     18 |        0       18 |     0     0 |     0     0 |      0     2       |
 | -3.95574e+37  1.60464e+07 | -3.95574e+37  1.60464e+07 | -3.95574e+37  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     0     0 |      2     0       |
 |  2.26972e+07  6.42211e+07 |  2.26972e+07  6.42211e+07 |  1.60464e+07  6.42211e+07 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL    6
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT         2.26972e+07
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     1
 SHIFT LEFT          1.67803e+07
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT           2.26972e+07
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      1
 RIGHT BOUNDARY OF THE INTERVAL   2.26972e+07
 LEFT BOUNDARY OF THE INTERVAL    1.60464e+07
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |     18 |        0       18 |     0     0 |     0     0 |      0     2       |
 | -3.95574e+37  1.60464e+07 | -3.95574e+37  1.60464e+07 | -3.95574e+37  3.21046e+07 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       23       24 |     0     0 |     1     0 |      2     0       |
 |  2.26972e+07  6.42211e+07 |  2.26972e+07  6.42211e+07 |  1.60464e+07  6.42211e+07 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    4
 ALL THE EIGENVALUES REQUESTED ON THIS LANCZOS RUN ARE FOUND
  
 NUMBER OF ACCEPTED EIGENVALUES   24
  
 ACCEPTED EIGENVALUES
  -1.16568e-09 -6.58454e-10 -3.27990e-10 -6.41666e-11  8.25525e-10  1.37336e-09
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
 |  +  |     24 |        0       24 |     0     0 |     0     0 |      0     0       |
 | -3.95574e+37  6.42211e+07 | -3.95574e+37  6.42211e+07 | -3.95574e+37  6.42211e+07 |
 -------------------------------------------------------------------------------------
  
 FIRST EIGENVALUE REQUESTED BY THE USER     1
  
 EIGENVALUES REQUESTED BY THE USER
  -1.16568e-09 -6.58454e-10 -3.27990e-10 -6.41666e-11  8.25525e-10  1.37336e-09
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04  9.29783e+04
   1.97239e+07  1.97239e+07  1.97239e+07  1.97239e+07  1.97239e+07  6.41026e+07
  
 EIGENVALUE EXTRACTION USING THE LANCZOS EIGENSOLVER IS COMPLETE
 NUMBER OF LANCZOS RUNS       :  42
 NUMBER OF LANCZOS STEPS      :  85
 NUMBER OF FACTORIZATIONS     :  34
 NUMBER OF BACK SUBSTITUTIONS :  96

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.34000    
       ELAPSED SYSTEM TIME (SEC)    =  8.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.42000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
1

   Abaqus 3DEXPERIENCE R2019x                                  Date 18-Mar-2020   Time 10:21:01
   For use by WILLIAM MARSH RICE UNIVERSITY INC under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    2     INCREMENT     1     STEP TIME    0.00    
 


                        S T E P       2     M A T R I X   G E N E R A T E
 
  THE FOLLOWING GLOBAL MATRICES WILL BE GENERATED AND WRITTEN TO Job-1_X2.sim

    STIFFNESS
    MASS
 
  THE MASS MATRIX PROPERTIES: CONSISTENT MASS MATRIX
 
  IF MULTY-POINT CONSTRAINTS ARE SPECIFIED IN THE MODEL, THEY WILL BE APPLIED
 
  THE FOLLOWING GLOBAL MATRICES WILL BE WRITTEN TO THE FOLLOWING TEXT FILES

    STIFFNESS          TO Job-1_STIF2.mtx
    MASS               TO Job-1_MASS2.mtx
 
  MATRIX OUTPUT FORMAT - MATRIX INPUT


          THE ANALYSIS HAS BEEN COMPLETED



     ANALYSIS SUMMARY:
     TOTAL OF          2  INCREMENTS
                       0  CUTBACKS IN AUTOMATIC INCREMENTATION
                       0  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                     130  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                      34  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       2  WARNING MESSAGES DURING USER INPUT PROCESSING
                       1  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                       0  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       0  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =   1.1000    
       SYSTEM TIME (SEC)    =  0.13000    
       TOTAL CPU TIME (SEC) =   1.2300    
       WALLCLOCK TIME (SEC) =          1
