# Simu – Simulation de dynamique de la chromatine et structures de type "beads"

Ce dépôt regroupe deux versions du code de simulation :
- **Simu_ring** : modélisation des *anneaux de chromatine* avec RNAP.
- **Simu_bead** : version adaptée pour des *grosses billes* interagissant par forces de Lennard-Jones et ressorts harmoniques.

---

## 📁 Structure du dépôt

Simu/
├── Simu_ring/
│   └── SImu_RNAP_chromatine/
│        └── N_1000/
│            ├── main_RNAP.c
│            ├── simulation.c
│            └── …
│
├── Simu_bead/
│   └── Bead_10/
│        └── N_1000/
│            ├── main_RNAP.c
│            ├── potentiels.c
│            └── …
│
└── .gitignore
