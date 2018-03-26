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
DIR="tests/Trpcage_amber_3"
cd ../../
nohup python main.py -input_file $DIR/ga_input.in &> $DIR/ga_run_script.out&
#python main.py -input_file  $DIR/ga_input.in > $DIR/ga_run_script.out
#python main.py -input_file $DIR/ga_input.in 2>&1 | tee $DIR/ga_run_script.out
#python main.py -input_file $DIR/ga_input.in

echo $! > $DIR/save_pid.txt
tail -f $DIR/ga_run_script.out
