#!/bin/bash
#PBS -N lemon_tight_clustMDD
#PBS -l nodes=1:ppn=1
#PBS -l walltime=4:00:00
#PBS -m abe

cd $PBS_O_WORKDIR
module load Java/1.8.0_231

#task tight clusters 
DIRECTORY_FILE="/scratch/gent/vo/000/gvo00027/projects/CBIGR/software/lemontree/"
java -jar ${DIRECTORY_FILE}lemontree_v3.1.1.jar -task tight_clusters -data_file Processed_counts_MDD.txt -cluster_file cluster_file_listMDD.txt -output_file LemonResults/tight_clustersMDD.txt -node_clustering true

echo End job 
