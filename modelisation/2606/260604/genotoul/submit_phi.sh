#!/bin/bash
#SBATCH -J phi_optim
#SBATCH -p workq
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH -o /work/user/ldodin/modelo/2606/260604/logs/job_%A_%a.out
#SBATCH -e /work/user/ldodin/modelo/2606/260604/logs/job_%A_%a.err

module purge
module load devel/Julia/1.11.4

WORKDIR=/work/user/ldodin/modelo/2606/260604
LOGDIR=$WORKDIR/logs
OUTDIR=$WORKDIR/output

mkdir -p "$LOGDIR" "$OUTDIR"

TASK_OFFSET=${TASK_OFFSET:-0}
REAL_ID=$(( SLURM_ARRAY_TASK_ID + TASK_OFFSET ))

COMBO=$(sed -n "${REAL_ID}p" $WORKDIR/../260603/combinations_15.txt)

echo "=========================================="
echo "Task   : $SLURM_ARRAY_TASK_ID / $SLURM_ARRAY_TASK_MAX"
echo "Node   : $SLURMD_NODENAME"
echo "Combo  : $REAL_ID / 3125"
echo "N_INTERIOR : $COMBO"
echo "Start  : $(date)"
echo "=========================================="

julia --project=/work/user/ldodin/modelo \
      --threads=$SLURM_CPUS_PER_TASK \
      $WORKDIR/fit_phi.jl $COMBO

echo "Done : $(date)"
