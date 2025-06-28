# Quantum ESPRESSO Workflow for Band Structure and DOS

This repository contains scripts to automate **Density Functional Theory (DFT)** calculations using [Quantum ESPRESSO](https://www.quantum-espresso.org/), focusing on electronic **band structure** and **density of states (DOS)** analysis.

---

## 📦 Contents

- `BANDS.sh` – Automates the band structure calculation.
- `DOS.sh` – Automates the total and projected DOS calculation.
- `plot_bands_dos.py` – Python script for plotting results using orbital projections.

---

## ✅ Requirements

- Quantum ESPRESSO ≥ 6.x  
- Bash shell  
- Python 3 with:
  - `numpy`
  - `matplotlib`

---

## ⚙️ Step-by-Step Instructions

### 1️⃣ Band Structure: `BANDS.sh`

1. Create a folder named `BANDS/` and place the necessary input files.
2. Run the band structure script:
   ```bash
   ./BANDS.sh
After running projwfc.x, open the file NbMoCO2proj.out.projwfc_up to inspect the projected orbital contributions.

For example, for Nb you may find:

4S 4P 4P 4P 5S 4D 4D 4D 4D 4D
0  1  2  3  4  5  6  7  8  9
Do the same for all atoms to identify the orbital indices. Example:

[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],         # Nb
[10, 11, 12, 13, 14, 15, 16, 17, 18, 19], # Mo
[20, 21, 22, 23],                        # C
[24, 25, 26, 27, 29, 30, 31, 28]         # O
These indices are used later for orbital-projected band structure plotting.

### 2️⃣ Density of States: `DOS.sh` 
Create a separate folder named DOS/ with input files.

Run the DOS script:

bash DOS.sh


### 3️⃣ Plotting
Gather the following files into a single folder:

atom_Mo.dat
atom_Nb.dat
atom_C.dat
atom_O.dat
NbMoCO2.k.pdos_tot
NbMoCO2proj.out.projwfc_up

Run the Python script:

python3 plot_bands_dos.py
