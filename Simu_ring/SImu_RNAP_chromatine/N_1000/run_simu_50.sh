#!/bin/bash
#SBATCH --job-name=simu_R1_n20_full
#SBATCH --output=slurm_simu_rnap_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=25          # Nombre de cœurs OpenMP disponibles
#SBATCH --mem=5G
#SBATCH -t 48:00:00
#SBATCH -p amd32

# Aller dans le dossier depuis lequel le job a été soumis
cd "$SLURM_SUBMIT_DIR"

module load openmpi4/4.1.1      # <- optionnel si pas d'MPI

# Compilation du code
echo "📦 Compilation du code..."
# gcc -g -O3 -ffast-math \
#     *.c -o main1 \
#     -lm -fopenmp \
#     -Wall \
#     -Wno-unused-variable \
#     -Wno-unused-but-set-variable \
#     -Wno-return-type \
#     -Wno-maybe-uninitialized \
#     -Wno-unused-result \
#     -Wno-comment

gcc- -g -O3 -ffast-math \
    main_RNAP.c \
    basic_functions.c \
    simulation.c \
    config.c \
    file.c \
    movement.c \
    neighborlist.c \
    potentiels.c \
    structures_depart.c \
    transcription_erwan.c \
    /Users/erwan/Documents/These/MTwister/mt19937ar.c \
    -Iinclude -lm -o \
    main \
    -lm -fopenmp \
    -Wall \
    -Wno-unused-variable \
    -Wno-unused-but-set-variable \
    -Wno-return-type \
    -Wno-maybe-uninitialized \
    -Wno-unused-result \
    -Wno-comment

echo "✅ Compilation OK"
chmod +x main1

# Threads OpenMP = cœurs SLURM
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

# Détection du parallélisme de fond (nombre de jobs simultanés lancés en &)
if [[ -n "${SLURM_CPUS_PER_TASK:-}" && "${SLURM_CPUS_PER_TASK}" -gt 1 ]]; then
  MAX_PARALLEL="${SLURM_CPUS_PER_TASK}"
else
  # Fallback : 1 si info indisponible
  MAX_PARALLEL=1
fi

echo "🖥️  Cœurs SLURM (OMP_NUM_THREADS) : $OMP_NUM_THREADS"
echo "⚙️  Jobs en parallèle (bg)        : $MAX_PARALLEL"

# Chemin absolu de l'exécutable
main1_path="$SLURM_SUBMIT_DIR/main1"

# Répertoire de base pour les résultats
parent_folder="$SLURM_SUBMIT_DIR/Simulations"
mkdir -p "$parent_folder"

# ======================
# Paramètres de simulation
# ======================
nb_rnap_values=(50)
vitesse_rnap_values=(0.1)
Ktranspt_values=(2)
seeds=(1 2 3 4 5)

running=0

# ======================
# Boucles de simulations
# ======================
for nb_rnap in "${nb_rnap_values[@]}"; do
  nb_rnap_folder="$parent_folder/nb-rnap_${nb_rnap}"
  mkdir -p "$nb_rnap_folder"

  if [[ "$nb_rnap" -ne 0 ]]; then
    for vitesse_rnap in "${vitesse_rnap_values[@]}"; do
      for Ktranspt in "${Ktranspt_values[@]}"; do
        k_folder="$nb_rnap_folder/vitesse_${vitesse_rnap}/Ktranspt_${Ktranspt}"
        mkdir -p "$k_folder"

        for seed in "${seeds[@]}"; do
          seed_folder="$k_folder/simulation_seed_${seed}"
          mkdir -p "$seed_folder"

          echo "▶ Simulation : nb_rnap=$nb_rnap | vitesse=$vitesse_rnap | Ktranspt=$Ktranspt | seed=$seed"
          (
            cd "$seed_folder"
            # Passage des 4 arguments : nb_rnap vitesse Ktranspt seed
            "$main1_path" "$nb_rnap" "$vitesse_rnap" "$Ktranspt" "$seed" > output.txt 2> error.txt
          ) &

          ((running++))
          if (( running >= MAX_PARALLEL )); then
            wait
            running=0
          fi
        done
      done
    done
  else
    # Cas témoin nb_rnap = 0 : on force vitesse=0 et Ktranspt=0
    for seed in "${seeds[@]}"; do
      seed_folder="$nb_rnap_folder/simulation_seed_${seed}"
      mkdir -p "$seed_folder"

      echo "▶ Simulation : nb_rnap=0 | vitesse=0 | Ktranspt=0 | seed=$seed"
      (
        cd "$seed_folder"
        "$main1_path" 0 0 0 "$seed" > output.txt 2> error.txt
      ) &

      ((running++))
      if (( running >= MAX_PARALLEL )); then
        wait
        running=0
      fi
    done
  fi
done

# Attendre la fin de toutes les tâches en arrière-plan
wait
echo "✅ Toutes les simulations sont terminées."
