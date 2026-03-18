#!/bin/bash
#SBATCH --job-name=GBSA
#SBATCH --output=GBSA_%A_%a.out
#SBATCH --error=GBSA_%A_%a.err
#SBATCH --partition=gpunode
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1

module purge
module load amber/24

$AMBERHOME/bin/MMPBSA.py -O -i mmgbsa_2_igb2.in -o gbsa_info.dat -sp solucion.prmtop -cp soluto.prmtop -rp complejo.prmtop -lp ligando.prmtop -y md.traj
