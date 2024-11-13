#!/bin/sh 
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J orig20_g0
### -- ask for number of cores (default: 1) -- 
#BSUB -n 16
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 18:00 
### -- specify that we need 2GB of memory per core/slot -- 
#BSUB -R "rusage[mem=100MB]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot -- 
#BSUB -M 3GB
### -- set the email address -- 
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
#BSUB -u s127782@student.dtu.dk
### -- send notification at start -- 
#BSUB -B 
### -- send notification at completion-- 
#BSUB -N 
### -- Specify the output and error file. %J is the job-id -- 
### -- -o and -e mean append, -oo and -eo mean overwrite -- 
#BSUB -o Output_%J.out 
#BSUB -e Error_%J.err 

# -- load the OpenFOAM module --
module load OpenFoam/v1912/gcc-9.2.0-openmpi-3.1.5

# -- program invocation here -- 
blockMesh
decomposePar
mpirun -np 16 DPMFoam -parallel > Output.log
