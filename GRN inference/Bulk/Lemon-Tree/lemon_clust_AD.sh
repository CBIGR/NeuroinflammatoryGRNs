#!/bin/bash
#PBS -N lemon_clust_AD
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -m abe

cd $PBS_O_WORKDIR
# call program
module load Java/1.8.0_231
#java -jar /scratch/gent/vo/000/gvo00027/projects/CBIGR/software/lemontree/lemontree_v3.1.1.jar
DIRECTORY_FILE="/scratch/gent/vo/000/gvo00027/projects/CBIGR/software/lemontree/"

# task ganesh
java -jar ${DIRECTORY_FILE}lemontree_v3.1.1.jar -task ganesh -data_file $3 -output_file ./LemonResults/$5

echo End job
