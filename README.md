# Quantum ESPRESSO Workflow for Band Structure and DOS

This repository contains scripts to automate **Density Functional Theory (DFT)** calculations using [Quantum ESPRESSO](https://www.quantum-espresso.org/), focusing on electronic **band structure** and **density of states (DOS)** analysis.

---

## Contents

- `BANDS.sh` – Automates the band structure calculation.
- `DOS.sh` – Automates the total and projected DOS calculation.
- `plot_bands_dos.py` – Python script for plotting results using orbital projections.

---

### Band Structure: `BANDS.sh`

1. Create a folder named `BANDS/` and place the necessary input files for the band structure calculation.
2. Run the band structure script:
   ```bash
   bash BANDS.sh
```

3. After running `projwfc.x`, open the file `NbMoCO2proj.out.projwfc_up` to inspect the projected orbital contributions.

   For example, for Nb, you may find:

   ```mathematica
   4S 4P 4P 4P 5S 4D 4D 4D 4D 4D
   0  1  2  3  4  5  6  7  8  9
   ```

   This indicates that the orbital indices for Nb are:

   ```python
   [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],         # Nb
   [10, 11, 12, 13, 14, 15, 16, 17, 18, 19], # Mo
   [20, 21, 22, 23],                        # C
   [24, 25, 26, 27, 29, 30, 31, 28]         # O
   ```

   These indices are used later for orbital-projected band structure plotting.

### Density of States: `DOS.sh`

1. Create a separate folder named `DOS/` and place the necessary input files for the DOS calculation.
2. Run the DOS script:

   ```bash
   bash DOS.sh
   ```
3. Once finished, extract the atomic-projected DOS using:

   ```bash
   sumpdos.x *\(Mo\)* > atom_Mo.dat
   sumpdos.x *\(Nb\)* > atom_Nb.dat
   sumpdos.x *\(C\)* > atom_C.dat
   sumpdos.x *\(O\)* > atom_O.dat
   ```

### Plotting

1. Gather the following files into a single folder:

   * `atom_Mo.dat`
   * `atom_Nb.dat`
   * `atom_C.dat`
   * `atom_O.dat`
   * `NbMoCO2.k.pdos_tot`
   * `NbMoCO2proj.out.projwfc_up`

2. Run the Python script to generate the plots:

   ```bash
   python3 plot_bands_dos.py
   ```

