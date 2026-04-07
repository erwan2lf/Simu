#!/usr/bin/env bash
set -euo pipefail

# ======================================================================
# CONFIG
# ======================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Détection des cores
if command -v nproc &> /dev/null; then
  NCORES=$(nproc)
elif [[ "$(uname)" == "Darwin" ]]; then
  NCORES=$(sysctl -n hw.ncpu)
else
  NCORES=4
fi

MAX_PARALLEL=$(( NCORES > 1 ? NCORES - 1 : 1 ))
echo "🖥️  Cœurs détectés      : $NCORES"
echo "⚙️  Jobs en parallèle    : $MAX_PARALLEL"

# ======================================================================
# COMPILATION
# ======================================================================
echo "📦 Compilation du code…"

gcc-14 -g -O3 -ffast-math -DMAC \
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
echo ""

# ======================================================================
# PARAMÈTRES DES SIMULATIONS
# ======================================================================

RESULT_DIR="$SCRIPT_DIR/Simulations"
mkdir -p "$RESULT_DIR"

nb_rnap_values=(10)
vitesse_rnap_values=(0.1)
Ktranspt_values=(2)
seeds=({1..2})

# ======================================================================
# FONCTION : lancer une seule fois la simulation
# ======================================================================
run_one_simulation() {
  local sim_dir="$1"
  local nb_rnap="$2"
  local vitesse="$3"
  local Ktranspt="$4"
  local seed="$5"

  cd "$sim_dir"

  echo "▶ Lancement simulation : nb_rnap=$nb_rnap | vitesse=$vitesse | Ktranspt=$Ktranspt | seed=$seed"

  echo ">> $(date +"[%H:%M:%S]") Lancement main (nb_rnap=$nb_rnap vitesse=$vitesse Ktranspt=$Ktranspt seed=$seed)" \
    | tee -a output.txt

  "$SCRIPT_DIR/main" "$nb_rnap" "$vitesse" "$Ktranspt" "$seed" >> output.txt 2>> error.txt

  echo "🏁 Fin du run" | tee -a output.txt

  cd "$SCRIPT_DIR"
}

# ======================================================================
# LANCEMENT DES SIMULATIONS
# ======================================================================

running=0

for nb_rnap in "${nb_rnap_values[@]}"; do

  nb_dir="$RESULT_DIR/nb-rnap_${nb_rnap}"
  mkdir -p "$nb_dir"

  for vitesse in "${vitesse_rnap_values[@]}"; do
    for Ktranspt in "${Ktranspt_values[@]}"; do

      if [[ "$nb_rnap" -eq 0 ]]; then
          kt_dir="$nb_dir/vitesse_${vitesse}"
          mkdir -p "$kt_dir"
      else
          kt_dir="$nb_dir/vitesse_${vitesse}/Ktranspt_${Ktranspt}"
          mkdir -p "$kt_dir"
      fi

      for seed in "${seeds[@]}"; do
          sim_dir="$kt_dir/simulation_seed_${seed}"
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
echo "🎉 Toutes les simulations sont terminées."