#!/bin/bash
#SBATCH --job-name=n0conf
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

source /opt/ohpc/admin/lmod/lmod/init/bash
module load openmpi4/4.1.1
module load fftw/3.3.8


# =============================
# Compilation du code
# =============================
echo "📦 Compilation du code..."

rm -f main *.o

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
    traj_binary.c \
    msd.c \
    forces.c \
    -Iinclude \
    -I${FFTW_DIR}/include \
    -L${FFTW_DIR}/lib \
    -Wl,-rpath,${FFTW_DIR}/lib \
    -lfftw3 -lm -fopenmp -o main \
    -Wall \
    -Wno-unused-variable \
    -Wno-unused-but-set-variable \
    -Wno-return-type \
    -Wno-maybe-uninitialized \
    -Wno-unused-result \
    -Wno-comment

if [ $? -ne 0 ]; then
    echo "❌ Compilation échouée — arrêt du job"
    exit 1
fi

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

nb_rnap_values=(0)
vitesse_rnap_values=(0.1)
Ktranspt_values=(2)
Delta_values=(1e-4)
gamma_fric_values=(1)
r_conf_values=(1000 15.88 13.6 11.6)
nb_passages_values=(1 2 4)         # ← nouveau
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
  local Delta="$6"
  local gamma_fric="$7"
  local r_conf="$8"
  local nb_passages="$9"           # ← nouveau

  cd "$sim_dir" || exit 1

  echo "▶ Simulation : nb_rnap=$nb_rnap | vitesse=$vitesse | K=$Ktranspt | seed=$seed | Delta=$Delta | gamma=$gamma_fric | r_conf=$r_conf | nb_passages=$nb_passages"
  echo ">> $(date +"[%H:%M:%S]") Lancement main" | tee -a output.txt

  "$SLURM_SUBMIT_DIR/main" "$nb_rnap" "$vitesse" "$Ktranspt" "$seed" "$Delta" "$gamma_fric" "$r_conf" "$nb_passages" \
      >> output.txt 2>> error.txt

  # --- Vérification trajectoire.bin ---
  if [ -f "Resultats/trajectoire.bin" ]; then
      echo "✅ trajectoire.bin présent ($(du -sh Resultats/trajectoire.bin | cut -f1))" \
          | tee -a output.txt
  else
      echo "⚠️  trajectoire.bin absent" | tee -a output.txt
  fi

  # --- Vérification msd.txt ---
  if [ -f "Resultats/msd.txt" ]; then
      n_lines=$(grep -c "^[^#]" Resultats/msd.txt 2>/dev/null || echo "?")
      echo "✅ msd.txt présent ($n_lines lags)" | tee -a output.txt
  else
      echo "⚠️  msd.txt absent" | tee -a output.txt
  fi

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
      for Delta in "${Delta_values[@]}"; do
        for gamma_fric in "${gamma_fric_values[@]}"; do
          for r_conf in "${r_conf_values[@]}"; do
            for nb_passages in "${nb_passages_values[@]}"; do   # ← nouveau

              k_folder="$nb_rnap_folder/vitesse_${vitesse}/Ktranspt_${Ktranspt}/Delta_${Delta}/gamma_${gamma_fric}/r_conf_${r_conf}/nb-passages_${nb_passages}"
              mkdir -p "$k_folder"

              for seed in "${seeds[@]}"; do
                sim_dir="$k_folder/simulation_seed_${seed}"
                mkdir -p "$sim_dir"

                (
                  run_one_simulation "$sim_dir" "$nb_rnap" "$vitesse" "$Ktranspt" "$seed" "$Delta" "$gamma_fric" "$r_conf" "$nb_passages"
                ) &

                ((running++))
                if (( running >= MAX_PARALLEL )); then
                  wait
                  running=0
                fi
              done

            done   # fin nb_passages
          done     # fin r_conf
        done
      done
    done
  done
done

wait

echo "======== END OF JOB $(date) ========"