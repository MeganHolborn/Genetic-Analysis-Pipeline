#!/usr/bin/env bash
#PBS -q long
#PBS -l walltime=900:00:00
#PBS -l nodes=1:ppn=1
#PBS -k oe
#PBS -N Snakemake

module load python-3.8.2
cd "$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
snakemake --cluster-config config/cluster.json --profile config/PBS-Torque-Profile 