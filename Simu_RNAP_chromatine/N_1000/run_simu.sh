#!/bin/bash
#SBATCH --job-name=c9.86_10r
#SBATCH --output=slurm_simu_rnap_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=5G
#SBATCH -t 48:00:00
#SBATCH -p amd32

echo "======== SLURM SIMULATION LAUNCH ========"
echo "Job ID: $SLURM_JOB_ID"
echo "Date  : $(date)"
echo "Dir   : $SLURM_SUBMIT_DIR"

cd "$SLURM_SUBMIT_DIR"

module load openmpi4/4.1.1

# =============================
# Compilation du code
# =============================
echo "📦 Compilation du code..."

gcc -g -O3 -ffast-math -DCLUSTER \
    main.c \
    basic_functions.c \
    simulation.c \
    config.c \
    file.c \
    movement.c \
    neighborlist.c \
    potentiels.c \
    structures_depart.c \
    transcription.c \
    mt19937ar.c \
    -Iinclude -lm -fopenmp -o main \
    -Wall \
    -Wno-unused-variable \
    -Wno-unused-but-set-variable \
    -Wno-return-type \
    -Wno-maybe-uninitialized \
    -Wno-unused-result \
    -Wno-comment

echo "✅ Compilation OK"
chmod +x main

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
MAX_PARALLEL=50

echo "🖥️  Threads OpenMP         : $OMP_NUM_THREADS"
echo "⚙️  Jobs en parallèle (bg)  : $MAX_PARALLEL"

# =============================
# PARAMÈTRES DES SIMULATIONS
# =============================
parent_folder="$SLURM_SUBMIT_DIR/Simulations"
mkdir -p "$parent_folder"

nb_rnap_values=(10)
vitesse_rnap_values=(0.1)
Ktranspt_values=(2)
seeds=({1..50})

# =============================
# Fonction : une simulation
# =============================
run_one_simulation() {
  local sim_dir="$1"
  local nb_rnap="$2"
  local vitesse="$3"
  local Ktranspt="$4"
  local seed="$5"

  cd "$sim_dir" || exit 1

  echo "▶ Simulation : nb_rnap=$nb_rnap | vitesse=$vitesse | K=$Ktranspt | seed=$seed"

  echo ">> $(date +"[%H:%M:%S]") Lancement main" | tee -a output.txt

  "$SLURM_SUBMIT_DIR/main" "$nb_rnap" "$vitesse" "$Ktranspt" "$seed" \
      >> output.txt 2>> error.txt

  echo "🏁 Fin du run (quel que soit le statut)" | tee -a output.txt

  cd "$SLURM_SUBMIT_DIR" || exit 1
}

# =============================
# Boucles des simulations
# =============================
running=0

for nb_rnap in "${nb_rnap_values[@]}"; do
  nb_rnap_folder="$parent_folder/nb-rnap_${nb_rnap}"
  mkdir -p "$nb_rnap_folder"

  for vitesse in "${vitesse_rnap_values[@]}"; do
    for Ktranspt in "${Ktranspt_values[@]}"; do
      k_folder="$nb_rnap_folder/vitesse_${vitesse}/Ktranspt_${Ktranspt}"
      mkdir -p "$k_folder"

      for seed in "${seeds[@]}"; do
        sim_dir="$k_folder/simulation_seed_${seed}"
        mkdir -p "$sim_dir"

        (
          run_one_simulation "$sim_dir" "$nb_rnap" "$vitesse" "$Ktranspt" "$seed"
        ) &

        ((running++))
        if (( running >= MAX_PARALLEL )); then
          wait
          running=0
        fi
      done
    done
  done
done

wait

echo "======== END OF JOB $(date) ========"
