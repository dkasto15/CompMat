#!/usr/bin/env bash
#SBATCH -p hebbe
#SBATCH -A SNIC2017-1-649
#SBATCH -J jobname
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -t 00:10:00
#SBATCH -o stdout
#SBATCH -e stderr

module purge
module load intel/2017b GPAW/1.3.0-Python-2.7.14

export GPAW_SETUP_PATH=$GPAW_SETUP_PATH:/c3se/apps/Glenn/gpaw/gpaw-setups-0.9.11271/
export GPAW_SETUP_PATH=$GPAW_SETUP_PATH:./

mpirun gpaw-python ./1_converge_kpoints_bulk0.py
