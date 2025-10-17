#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Détection du nombre de cœurs
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

# (optionnel) limiter OpenMP
# export OMP_NUM_THREADS=1

# Compilation
echo "📦 Compilation du code…"
gcc-14 -std=c11 -g -O3 -ffast-math \
  *.c -o "$SCRIPT_DIR/main1" -lm -fopenmp \
  -Wall \
  -Wno-unused-variable \
  -Wno-unused-but-set-variable \
  -Wno-return-type \
  -Wno-maybe-uninitialized \
  -Wno-unused-result \
  -Wno-comment
echo "✅ Compilation OK"

# Paramètres
RESULT_DIR="$SCRIPT_DIR/Simulations"
mkdir -p "$RESULT_DIR"

nb_rnap_values=(1 3 5)
vitesse_rnap_values=(0.1)
Ktranspt_values=(8 16 64)
alpha_values=(2 4 8)         
seeds=(1 2 3 4 5)

running=0

# Boucle de simulations
for nb_rnap in "${nb_rnap_values[@]}"; do
  nb_dir="$RESULT_DIR/nb-rnap_${nb_rnap}"
  mkdir -p "$nb_dir"

  if (( nb_rnap > 0 )); then
    for vitesse in "${vitesse_rnap_values[@]}"; do
      for Ktranspt in "${Ktranspt_values[@]}"; do
        for alpha in "${alpha_values[@]}"; do
          combo_dir="$nb_dir/vitesse_${vitesse}/Ktranspt_${Ktranspt}/alpha_${alpha}"
          mkdir -p "$combo_dir"

          for seed in "${seeds[@]}"; do
            sim_dir="$combo_dir/simulation_seed_${seed}"
            mkdir -p "$sim_dir"

            echo "▶ Lancement : nb_rnap=$nb_rnap | vitesse=$vitesse | Ktranspt=$Ktranspt | alpha=$alpha | seed=$seed"
            (
              cd "$sim_dir"
              # main1 args: nb_rnap vitesse Ktranspt alpha seed
              "$SCRIPT_DIR/main1" "$nb_rnap" "$vitesse" "$Ktranspt" "$alpha" "$seed" \
                >> "output.txt" 2>> "error.txt"
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
  else
    # Cas sans RNAP : on peut fixer vitesse=0 Ktranspt=0 alpha=0
    for seed in "${seeds[@]}"; do
      sim_dir="$nb_dir/vitesse_0/Ktranspt_0/alpha_0/simulation_seed_${seed}"
      mkdir -p "$sim_dir"

      echo "▶ Lancement : nb_rnap=0 | vitesse=0 | Ktranspt=0 | alpha=0 | seed=$seed"
      (
        cd "$sim_dir"
        "$SCRIPT_DIR/main1" 0 0 0 0 "$seed" \
          >> "output.txt" 2>> "error.txt"
      ) &

      ((running++))
      if (( running >= MAX_PARALLEL )); then
        wait
        running=0
      fi
    done
  fi
done

wait
echo "✅ Toutes les simulations sont terminées."