kmeans clustering for studying migrations in proteins
=====================================================

C++ implementation of the kmeans algorithm (http://en.wikipedia.org/wiki/K-means_clustering)
for finding clusters of stable microstates when studying ligand migration in proteins
 
Molecular systems are realigned using the Kabsch algorithm (http://en.wikipedia.org/wiki/Kabsch_algorithm)
 
Copyright (c) 2015, Florent Hédin, Pierre-André Cazade, and the University of Basel.
All rights reserved.
 
The 3-clause BSD license is applied to this software.
See LICENSE.txt

-----------
Compiling :
-----------

The Kabsch algorithm implementation for realigning the dcd frames requires the BLAS and LAPACK development packages to be installed on your system : 

On rpm based systems:

    sudo yum install blas-devel lapack-devel

On deb based systems:

    sudo apt-get install libblas-dev liblapack-dev

------------------
Print help/usage :
------------------

./kmeans

Error, not enough arguments, usage is : kmeans -idx {index of atom to study (taken from a PSF for example)} -dcd {number of dcd files} {paths to DCDs} -out {path to outputFile} -xyz {path to outputXYZ}
inputDCDs and outputFile and outputXYZ are necessary fileNames

optional arguments : 
-align {first} {last}    if present, align all frames from all dcds relative to first frame of first dcd before starting clustering
-interactive     if present the user will have to provide parameters interactively
-cycles          provide number of cycles
-cutoff          provide cutoff value for cluster determination
-mult    provide the Multiplicator factor of the cutoff for exclusion of lonely microstates
-thresh          provide the Threshold to consider a microstate is in water
-tolerance       provide the Tolerance for convergence

----------------------------------------------------
Run with all parameters provided from command line : 
----------------------------------------------------

If we want to use all coordinates of atom 2534 from 2 dcds as starting points for the list of mictro states, aligning all dcd frames to atoms 1 to 2532

./kmeans -idx 2536 -dcd 2 ../1111/run_0.dcd ../rst1/1111/rst1_0.dcd -out test1111.out -xyz test1111.xyz -cutoff 1.8 -mult 4.5 -align 1 2532



