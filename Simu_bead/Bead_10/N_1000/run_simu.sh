#!/bin/bash
#SBATCH --job-name=simu_R1_n20_full
#SBATCH --output=slurm_simu_rnap_%j.out
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=25                # Nombre de cœurs alloués au job
#SBATCH --mem=5G
#SBATCH -t 48:00:00
#SBATCH -p amd32

# Aller dans le dossier depuis lequel le job a été soumis
cd $SLURM_SUBMIT_DIR

# Chargement des modules nécessaires
module load openmpi4/4.1.1                # Adapte si nécessaire

# Compilation du code
echo "📦 Compilation du code..."
gcc -std=c11 -g -O3 -ffast-math *.c -o main1 -lm -fopenmp \
    -Wall \
    -Wno-unused-variable \
    -Wno-unused-but-set-variable \
    -Wno-return-type \
    -Wno-maybe-uninitialized \
    -Wno-unused-result \
    -Wno-comment

if [ $? -ne 0 ]; then
    echo "❌ Erreur de compilation"
    exit 1
fi

chmod +x main1

# Chemin absolu de l'exécutable
main1_path="$SLURM_SUBMIT_DIR/main1"

# Répertoire de base pour les résultats
parent_folder="$SLURM_SUBMIT_DIR/Simulations"
mkdir -p "$parent_folder"

# Paramètres de simulation
nb_rnap_values=(10)
vitesse_rnap_values=(0.1)
Ktranspt_values=(8 16)
alpha_values=(4)
seeds=(1 2 3 4 5)

# Nombre de tâches en parallèle autorisées
MAX_PARALLEL=$SLURM_CPUS_PER_TASK
running=0

# Boucle sur tous les paramètres
for nb_rnap in "${nb_rnap_values[@]}"; do
    nb_rnap_folder="$parent_folder/nb-rnap_$nb_rnap"
    mkdir -p "$nb_rnap_folder"

    if [ "$nb_rnap" -ne 0 ]; then
        for vitesse_rnap in "${vitesse_rnap_values[@]}"; do
            for Ktranspt in "${Ktranspt_values[@]}"; do
                for alpha in "${alpha_values[@]}"; do
                    combo_folder="$nb_rnap_folder/vitesse_${vitesse_rnap}/Ktranspt_${Ktranspt}/alpha_${alpha}"
                    mkdir -p "$combo_folder"

                    for seed in "${seeds[@]}"; do
                        seed_folder="$combo_folder/simulation_seed_$seed"
                        mkdir -p "$seed_folder"

                        echo "▶ Simulation : nb_rnap=$nb_rnap, vitesse=$vitesse_rnap, Ktranspt=$Ktranspt, alpha=$alpha, seed=$seed"
                        (
                            cd "$seed_folder" || exit 1
                            "$main1_path" "$nb_rnap" "$vitesse_rnap" "$Ktranspt" "$alpha" "$seed" \
                                > output.txt 2> error.txt
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
        for seed in "${seeds[@]}"; do
            seed_folder="$nb_rnap_folder/vitesse_0/Ktranspt_0/alpha_0/simulation_seed_$seed"
            mkdir -p "$seed_folder"

            echo "▶ Simulation : nb_rnap=0 (vitesse=0, Ktranspt=0, alpha=0), seed=$seed"
            (
                cd "$seed_folder" || exit 1
                "$main1_path" 0 0 0 0 "$seed" > output.txt 2> error.txt
            ) &

            ((running++))
            if (( running >= MAX_PARALLEL )); then
                wait
                running=0
            fi
        done
    fi
done

# Attendre que toutes les simulations soient terminées
wait
echo "✅ Toutes les simulations sont terminées."
