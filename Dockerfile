# EVOLVE Dockerfile for CI

# Github Actions does not support GPU acceleration and the parallized amber pmemd engine is not available for free via conda. 
# The unittests regarding this are not executed


FROM continuumio/miniconda3

# Setup a spot for the code
WORKDIR /evolve
# Install Python dependencies
COPY environment.yml environment.yml
RUN conda env create -f environment.yml

# Activate environment
SHELL ["conda", "run", "-n", "evolve", "/bin/bash", "-c"]


# Copy necessary files
COPY share share/
COPY unittests unittests/
COPY main.py .
COPY evaluators.py .

# This is necessary due to a bug in Docker
RUN true

COPY src ./src/

# Set Environment variables
ENV PYTHONPATH="/evolve/:${PYTHONPATH}"

ENV PATH="/opt/conda/envs/evolve/bin:$PATH"

ENV CONDA_DEFAULT_ENV evolve

ENV EVOLVE_GPU=FALSE

SHELL ["conda", "run", "-n", "evolve", "/bin/bash", "-c"]

# Chceck that openbabel is available
RUN echo "Make sure openbabel is installed:"
RUN python -c "from openbabel import openbabel"

# Change to unittest directory
WORKDIR /evolve/unittests


# Execute tests
RUN ["python", "test_gaapi.py"]
RUN ["python", "test_openbabel.py"]
RUN ["python", "test_StartGA_Run.py"]


RUN ["/bin/bash"]