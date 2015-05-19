kmeans clustering for studying migrations in proteins
=====================================================

C++ implementation of the kmeans algorithm (http://en.wikipedia.org/wiki/K-means_clustering)
for finding clusters of stable microstates when studying ligand migration in proteins
 
Copyright (c) 2015, Florent Hédin, Pierre-André Cazade, and the University of Basel.
All rights reserved.
 
The 3-clause BSD license is applied to this software.
See LICENSE.txt

Print help/usage :
------------------

./kmeans

Error, not enough arguments, usage is : ./kmeans -idx {index of atom to study (taken from a PSF for example)} -dcd {number of dcd files} {paths to DCDs} -out {path to outputFile} -xyz {path to outputXYZ}
inputFile and outputFile and outputXYZ are necessary fileNames

optional arguments : 
-interactive     if present the user will have to provide parameters interactively
-cycles          provide number of cycles
-cutoff          provide cutoff value for cluster determination
-mult    provide the Multiplicator factor of the cutoff for exclusion of lonely microstates
-thresh          provide the Threshold to consider a microstate is in water
-tolerance       provide the Tolerance for convergence



Run with all parameters provided from command line : 
----------------------------------------------------

If we want to use all coordinates of atom 2534 from 2 dcds as starting points for the list of mictro states : 

./kmeans -idx 2533 -dcd 2 ../../1111/run_0.dcd  ../1111/rst1_0.dcd -out outtest -xyz xyztest -cutoff 2.0 -mult 1.0

User provided 2 dcd files : 
../../1111/run_0.dcd
../1111/rst1_0.dcd
Values used for parameters :
         rCutoff : 2
         mult : 1
         rThrs : 25
         Tol : 0.0001
         maxCycle : 250
Global exclusion (Cutoff*MultiplicatorFactor) is rExclude = 2

Reading coordinates from dcd : ../../1111/run_0.dcd
HDR :   CORD
ICNTRL :        5000    10100   100     500000  0       0       0       36711   0       1017614562   0       0       0       0       0       0       0       0       0       40
NTITLE :        3
TITLE : * MBXE4                                                                         * PROJECT: RELAXATION TRAJECTORIES AFTER XE REMOVAL                             *  DATE:     5/ 7/15     17:10:34      CREATED BY USER: hedin                   
NATOM : 17725
LNFREAT :       17725
Done for dcd : ../../1111/run_0.dcd

Reading coordinates from dcd : ../1111/rst1_0.dcd
HDR :   CORD
ICNTRL :        5000    100     100     500000  0       0       0       36711   0       1017614562   0       0       0       0       0       0       0       0       0       40
NTITLE :        3
TITLE : * MBXE4                                                                         * PROJECT: RELAXATION TRAJECTORIES AFTER XE REMOVAL                             *  DATE:     5/12/15     17: 0:43      CREATED BY USER: hedin                   
NATOM : 17725
LNFREAT :       17725
Done for dcd : ../1111/rst1_0.dcd

Total number of frames read from the 2 dcds is : 10000

Now performing the clustering task ...

Not Converged! Number of clusters: 17 drMin: 0 drMax 1.00187
4 4250 42.5
9 1935 19.35
0 1875 18.75
11 759 7.59
12 405 4.05
8 284 2.84
3 81 0.81
14 70 0.7
10 58 0.58
16 51 0.51
7 49 0.49
13 46 0.46
1 37 0.37
6 34 0.34
2 31 0.31
5 23 0.23
15 12 0.12

