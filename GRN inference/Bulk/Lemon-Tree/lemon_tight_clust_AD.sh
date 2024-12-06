#!/bin/bash
#PBS -N lemon_tight_clustAD
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -m abe

cd $PBS_O_WORKDIR
module load Java/1.8.0_231

#task tight clusters 
DIRECTORY_FILE="/scratch/gent/vo/000/gvo00027/projects/CBIGR/software/lemontree/"
java -jar ${DIRECTORY_FILE}lemontree_v3.1.1.jar -task tight_clusters -data_file Processed_counts_AD.txt -cluster_file cluster_file_list.txt -output_file LemonResults/tight_clustersAD.txt -node_clustering true

echo End job 