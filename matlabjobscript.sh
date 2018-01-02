#!/bin/bash
#SBATCH -t 15:00:00
#SBATCH -n 10
#SBATCH -N 1
#SBATCH -A physics-hi
#SBATCH --mail-user=paraj@umd.edu
#SBATCH --mail-type=END
#SBATCH --mem-per-cpu=1024
#SBATCH --share
#SBATCH --array=1-1
. ~/.profile
module load matlab/2017a

echo "
clear

Delta=0.01; %define a time-step 
Tf=10^3; %Final time
T0=0;%Initial time
[x, w] = GetTrapRule((Tf-T0)/Delta+1, T0, Tf);
Ns = [100];% array of system sizes
Bs =2;% 1.8:0.05:2.2;
J = 1;
narr=floor(logspace(1,log10((Tf-T0)/Delta+1),25));%generate log-space array to get times.
ms=calculate_avgt(Ns, Bs, J,narr, x, w); 

JobID='$SLURM_JOB_ID';
">>TFIM_timedependence_params_$SLURM_JOB_ID.m

matlab -nosplash -nodesktop -nojvm < TFIM_timedependence_params$SLURM_JOB_ID.m > output_$SLURM_JOB_ID.txt
##rm *$SLURM_JOB_ID.m

hostname
date

