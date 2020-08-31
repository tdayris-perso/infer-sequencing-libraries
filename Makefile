SHELL := bash
.ONESHELL:
.SHELLFLAGS := -euic
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

### Variables ###
# Tools
PYTEST           = pytest
BASH             = bash
CONDA            = conda
PYTHON           = python3.8
SNAKEMAKE        = snakemake
CONDA_ACTIVATE   = source $$(conda info --base)/etc/profile.d/conda.sh && conda activate && conda activate

# Paths
TEST_PIPELINE    = scripts/prepare_pipeline.py
SNAKE_FILE       = Snakefile
ENV_YAML         = envs/workflow.yaml
READS_PATH       = '${PWD}/tests/reads'
FASTA_PATH       = '${PWD}/tests/genome/genome.chr21.fa'
GTF_PATH         = '${PWD}/tests/genome/annotation.chr21.gtf'
BED_PATH         = '${PWD}/tests/genome/chr21.bed'

# Arguments
ENV_NAME         = infer-sequencing-libraries
SNAKE_THREADS    = 1
PYTEST_ARGS      = -vv

# Recipes
default: all-unit-tests


all-unit-tests:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTEST} ${PYTEST_ARGS} ${TEST_PIPELINE}
.PHONY: all-unit-tests


# Environment building through conda
conda-tests:
	${CONDA_ACTIVATE} base && \
	${CONDA} env create --file ${ENV_YAML} --force && \
	${CONDA} activate ${ENV_NAME}
.PHONY: conda-tests


test-conda-report.html:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTHON} ${TEST_PIPELINE} ${READS_PATH} ${FASTA_PATH} ${GTF_PATH} ${BED_PATH} && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --printshellcmds --reason --forceall --configfile ${PWD}/config.yaml && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --printshellcmds --reason --forceall --configfile ${PWD}/config.yaml --report test-conda-report.html


clean:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --printshellcmds --reason --forceall --configfile ${PWD}/config.yaml --delete-all-output
.PHONY: clean



# Display pipeline graph
workflow.png:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTHON} ${TEST_PIPELINE} ${READS_PATH} ${FASTA_PATH} ${GTF_PATH} ${BED_PATH} && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --printshellcmds --reason --forceall --configfile ${PWD}/config.yaml --rulegraph | dot -T png > workflow.png

example.png:
	${CONDA_ACTIVATE} ${ENV_NAME} && \
	${PYTHON} ${TEST_PIPELINE} ${READS_PATH} ${FASTA_PATH} ${GTF_PATH} ${BED_PATH} && \
	${SNAKEMAKE} -s ${SNAKE_FILE} --use-conda -j ${SNAKE_THREADS} --printshellcmds --reason --forceall --configfile ${PWD}/config.yaml --dag | dot -T png > example.png
