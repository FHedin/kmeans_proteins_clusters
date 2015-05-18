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

Error, not enough arguments, usage is : kmeans -idx index of atom to study (taken from a PSF for example) -dcd {path to DCD} -out {path to outputFile} -xyz {path to outputXYZ}
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

If we want to use all coordinates of atom 2534 in dcd as tarting points for the list of mictro states : 

./kmeans -idx 2534 -dcd data_md/step_1.dcd -out out.dat -xyz out.xyz -cutoff 1.7 -mult 1.0 -thresh 25.0 -tolerance 0.0001

The corresponding PSF was :

* MBXE4
* PROJECT: RELAXATION TRAJECTORIES AFTER XE REMOVAL
...
...
    2533 XE3  1    XE3  XE   XE      0.00000       131.293           0   0.00000     -0.301140E-02
    2534 XE4  1    XE4  XE   XE      0.00000       131.293           0   0.00000     -0.301140E-02
...
...

i.e. we see that atoms 2533 and 2534 were XE atoms for which we want to perform kmeans analysis.

as output we get :

HDR :   CORD
ICNTRL :        1000    500     500     500000  0       0       0       36705   0       1017614562      0       0       0       0       0       0       0       0       0       40
NTITLE :        3
TITLE : * MBXE4                                                                         * PROJECT: RELAXATION TRAJECTORIES AFTER XE REMOVAL                             *  DATE:     4/13/15      9:53: 0      CREATED BY USER: hedin                   
NATOM : 17723
LNFREAT :       17723
Values used for parameters :
         rCutoff : 1.7
         mult : 1
         rThrs : 25
         Tol : 0.0001
         maxCycle : 250
Global exclusion (Cutoff*MultiplicatorFactor) is rExclude = 1.7
Not Converged! Number of clusters: 11 drMin: 0.124545 drMax 0.704248
1 262 26.2
0 198 19.8
7 139 13.9
2 97 9.7
9 97 9.7
3 75 7.5
4 62 6.2
6 33 3.3
8 24 2.4
5 7 0.7
10 6 0.6



the 6 first lines show some data read at the beginning of the dcd file.

Then come a summary of parameters used for clustering.
