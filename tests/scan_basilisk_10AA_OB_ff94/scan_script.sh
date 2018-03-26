#!/bin/bash
# bash script to run several times GA with or without using Basilisk

# load required modules to run GA
# module purge
source /software/ENV/modules
module load anaconda/4.3.1
module load openbabel/2.3.2/gcc-4.9.2
module load amber/16/gcc-6.3.0
module load g16


CURRENT_DIR="scan_basilisk_10AA_OB_ff94"
LIST="1 2 3 4 5 6 7 8 9 10"

# loop over simulations without basilisk
for i in $LIST
do

WORKING_DIR="no_basilisk_$i"
mkdir $WORKING_DIR
cd $WORKING_DIR

rm -f ga_input.in
cat > ga_input.in << EOF
# SAMPLE GABIO INPUT FILE 

[OPTIMIZATION]
composition_optimization = False
dihedral_optimization = True
dihedral_residue_indexes = 0 1 2 3 4 5 6 7 8 9
basilisk_and_sidechains = False

[MOLECULE]
initial_molecule = tests/$CURRENT_DIR/reduce_10modified.pdb

[IO]
output_path = tests/$CURRENT_DIR/$WORKING_DIR
output_file = OUT_molec.pdb
store_intermediate_mol = True

[GA_GLOBAL]
seed_population = False
population_file_path = None
population_size = 30 
max_iteration = 50

[EVALUATOR]
evaluators = openbabel_ff94

[MC_GENERATE]
mc_generate_dihedrals = True
distribution_path = DO_NOT_CHANGE
mcmove_lbound = -45.0
mcmove_ubound = 45.0
num_mc_steps = 1000
dihedral_probability_pointers = 0 0 0 0 0 0 0 0 1 0 

[GA_SELECTION]
selector = tournament_wor
tournament_size = 2

[GA_CROSSOVER]
crossover = sbx
mating_probability = 0.8
genewise_crossover_probability = 0.5
sbx_distribution_index = 5.0

[GA_MUTATION]
mutator = poly_dihedral
mutation_probability = 0.1
genewise__mutation_probability = 0.1
poly_eta = 5.0

[GA_REPLACER]
replacer = elitism
elitist_factor = 0.4

[DEBUG]
verbose = False
EOF

# Execute GA
cd ../../..
nohup python main.py -input_file tests/$CURRENT_DIR/$WORKING_DIR/ga_input.in &> tests/$CURRENT_DIR/$WORKING_DIR/ga_run_script.out &
echo $! > tests/$CURRENT_DIR/$WORKING_DIR/save_pid.txt
#tail -f tests/$CURRENT_DIR/$WORKING_DIR/ga_run_script.out
cd tests/$CURRENT_DIR

done


# loop over simulations with basilisk
for i in $LIST
do

WORKING_DIR="basilisk_$i"
mkdir $WORKING_DIR
cd $WORKING_DIR

rm -f ga_input.in
cat > ga_input.in << EOF
# SAMPLE GABIO INPUT FILE 

[OPTIMIZATION]
composition_optimization = False
dihedral_optimization = True
dihedral_residue_indexes = 0 1 2 3 4 5 6 7 8 9
basilisk_and_sidechains = True

[MOLECULE]
initial_molecule = tests/$CURRENT_DIR/reduce_10modified.pdb

[IO]
output_path = tests/$CURRENT_DIR/$WORKING_DIR
output_file = OUT_molec.pdb
store_intermediate_mol = True

[GA_GLOBAL]
seed_population = False
population_file_path = None
population_size = 30 
max_iteration = 50

[EVALUATOR]
evaluators = openbabel_ff94

[MC_GENERATE]
mc_generate_dihedrals = True
distribution_path = DO_NOT_CHANGE
mcmove_lbound = -45.0
mcmove_ubound = 45.0
num_mc_steps = 1000
dihedral_probability_pointers = 0 0 0 0 0 0 0 0 1 0 

[GA_SELECTION]
selector = tournament_wor
tournament_size = 2

[GA_CROSSOVER]
crossover = sbx
mating_probability = 0.8
genewise_crossover_probability = 0.5
sbx_distribution_index = 5.0

[GA_MUTATION]
mutator = poly_dihedral
mutation_probability = 0.1
genewise__mutation_probability = 0.1
poly_eta = 5.0

[GA_REPLACER]
replacer = elitism
elitist_factor = 0.4

[DEBUG]
verbose = False
EOF

# Execute GA
cd ../../..
nohup python main.py -input_file tests/$CURRENT_DIR/$WORKING_DIR/ga_input.in &> tests/$CURRENT_DIR/$WORKING_DIR/ga_run_script.out &
echo $! > tests/$CURRENT_DIR/$WORKING_DIR/save_pid.txt
#tail -f tests/$CURRENT_DIR/$WORKING_DIR/ga_run_script.out
cd tests/$CURRENT_DIR

done
