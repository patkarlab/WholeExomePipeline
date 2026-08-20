#!/bin/bash
#PBS -q short
#PBS -N lymphoma
#PBS -l select=2:ncpus=56:mem=250gb
#PBS -l walltime=24:00:00
#PBS -k od
#PBS -j oe
#PBS -o /home/patkarlab-clinical/pipelines/WholeExomePipeline/lymphoma.out

cd "/home/patkarlab-clinical/pipelines/WholeExomePipeline"

echo "========================================"
echo "Job ID: $PBS_JOBID"
echo "========================================"

# ---- Node info (PBS equivalent of SLURM)
NODE_NAME=$(hostname)
echo "Node: $NODE_NAME"

# Full node details
NODE_INFO=$(pbsnodes "$NODE_NAME")

# Total CPUs on node
CPUTOT=$(echo "$NODE_INFO" | awk -F= '/resources_available.ncpus/ {print $2}' | xargs)

# CPUs already allocated on node
CPUASSIGNED=$(echo "$NODE_INFO" | awk -F= '/resources_assigned.ncpus/ {print $2}' | xargs)
[ -z "$CPUASSIGNED" ] && CPUASSIGNED=0

echo "----------------------------------------"
echo "Total CPUs on node     : $CPUTOT"
echo "Allocated CPUs (node) : $CPUASSIGNED"
echo "----------------------------------------"


# ---- Run pipeline

./run_nextflow_lymphoma.sh > script.log
