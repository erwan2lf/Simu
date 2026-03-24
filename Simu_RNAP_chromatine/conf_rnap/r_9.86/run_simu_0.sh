#!/bin/bash
#SBATCH --job-name=c9_0r
#SBATCH --output=slurm_simu_rnap_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=5G
#SBATCH -t 48:00:00
#SBATCH -p amd32

cd "$SLURM_SUBMIT_DIR"

module load openmpi4/4.1.1

echo "📦 Compilation du code..."
echo "📂 Working directory = $(pwd)"
echo "📂 Contenu du dossier avant compilation :"
ls -l *.c include/*.h 2>/dev/null || echo "[DEBUG] Aucun fichier trouvé"

echo "🧾 Liste EXACTE des fichiers .c passés à gcc :"
cat <<EOF
main.c
basic_functions.c
simulation.c
config.c
file.c
movement.c
neighborlist.c
potentiels.c
structures_depart.c
transcription.c
mt19937ar.c
EOF

echo "🔎 Hash des fichiers source (pour vérifier la version utilisée) :"
md5sum \
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
    mt19937ar.c 2>/dev/null || echo "[DEBUG] md5sum indisponible"

echo "🚧 Lancement de gcc..."

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
echo "🔎 Infos sur le binaire généré :"
ls -l main
file main
strings main | grep -i cohes | head -n 10 || echo "[DEBUG] aucune trace de 'cohes' dans le binaire"

chmod +x main

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
MAX_PARALLEL=50

main_path="$SLURM_SUBMIT_DIR/main"
parent_folder="$SLURM_SUBMIT_DIR/Simulations"
mkdir -p "$parent_folder"

# ======================
# Paramètres de simulation
# ======================
nb_rnap_values=(0)
vitesse_rnap_values=(0.1)
Ktranspt_values=(4)
seeds=({1..50})

running=0

# ======================
# Boucles de simulations
# ======================
for nb_rnap in "${nb_rnap_values[@]}"; do
  nb_rnap_folder="$parent_folder/nb-rnap_${nb_rnap}"
  mkdir -p "$nb_rnap_folder"

  if [[ "$nb_rnap" -eq 0 ]]; then

    for vitesse_rnap in "${vitesse_rnap_values[@]}"; do
      vitesse_folder="$nb_rnap_folder/vitesse_${vitesse_rnap}"
      mkdir -p "$vitesse_folder"

      for seed in "${seeds[@]}"; do

        seed_folder="$vitesse_folder/simulation_seed_${seed}"
        mkdir -p "$seed_folder"

        echo "▶ Contrôle : nb_rnap=0 | vitesse=${vitesse_rnap} | seed=$seed"

        (
          cd "$seed_folder"
          "$main_path" 0 "$vitesse_rnap" 0 "$seed" > output.txt 2> error.txt
        ) &

        ((running++))
        if (( running >= MAX_PARALLEL )); then
          wait
          running=0
        fi
      done
    done

  else
    for vitesse_rnap in "${vitesse_rnap_values[@]}"; do
      for Ktranspt in "${Ktranspt_values[@]}"; do

        k_folder="$nb_rnap_folder/vitesse_${vitesse_rnap}/Ktranspt_${Ktranspt}"
        mkdir -p "$k_folder"

        for seed in "${seeds[@]}"; do

          seed_folder="$k_folder/simulation_seed_${seed}"
          mkdir -p "$seed_folder"

          echo "▶ Simulation : nb_rnap=$nb_rnap | vitesse=$vitesse_rnap | K=$Ktranspt | seed=$seed"

          (
            cd "$seed_folder"
            "$main_path" "$nb_rnap" "$vitesse_rnap" "$Ktranspt" "$seed" > output.txt 2> error.txt
          ) &

          ((running++))
          if (( running >= MAX_PARALLEL )); then
            wait
            running=0
          fi
        done
      done
    done
  fi
done

wait
echo "🎉 Toutes les simulations sont terminées."
