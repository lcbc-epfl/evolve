
source /software/ENV/modules

module load amber/18/intel-17.0.4
module load openbabel/2.4.1/gcc-6.3.0
module unload anaconda/5.0.1/python-3.6

module load cuda/10.1

source /home/duerr/miniconda3/bin/activate openmm
export PYTHONPATH=/home/duerr/phd/02_Projects/EVOLVE/Development/evolve/:$PYTHONPATH


python test_openbabel.py
python test_gaapi.py
python test_StartGA_Run.py
