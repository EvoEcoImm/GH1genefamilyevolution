#!/bin/bash

#SBATCH --job-name=iq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=2000
#SBATCH --time=6:00:00
#SBATCH --output=GH1s.genome.protein.final.faa.muscle.rascal.output.txt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=shulin.he@UCLouvain.be

iqtree -s GH1s.genome.protein.final.faa.muscle.rascal -m LG+R10 -nt 12 -bb 3000 -alrt 1000 -seed 7896788
