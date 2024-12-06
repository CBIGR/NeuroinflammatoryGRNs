#!/bin/bash
#PBS -N lemon_reg_MDD
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -m abe

cd $PBS_O_WORKDIR
module load Java/1.8.0_231

#task regulators 
DIRECTORY_FILE="/scratch/gent/vo/000/gvo00027/projects/CBIGR/software/lemontree/"
java -jar ${DIRECTORY_FILE}lemontree_v3.1.1.jar -task regulators -data_file Processed_counts_MDD.txt -reg_file RegulatorsMDD.txt -cluster_file LemonResults/tight_clustersMDD.txt -output_file LemonResults/reg_tfMDD

echo End job 
