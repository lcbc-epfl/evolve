#!/bin/bash
# Bash script to run GA

# Load required modules to run GA
# module purge
source /software/ENV/modules
module load anaconda/4.3.1
module load openbabel/2.3.2/gcc-4.9.2
module load amber/16/gcc-6.3.0
module load g16

# Execute GA
cd ../../
python main.py -input_file tests/amber_gpgg/ga_input.in 
