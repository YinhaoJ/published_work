# EInsights into Protein Unfolding under pH, Temperature, and Shear Using Molecular Dynamics Simulations
This repository contains input files, scripts, and protein structures used in the simulations for the following publication:

*Insights into Protein Unfolding under pH, Temperature, and Shear Using Molecular Dynamics Simulations*

---

## Content

### Protein Structure
- **7_4pH_BSA.pdb**  
  Structure of Bovine Serum Albumin at pH 7.4 (CHARMM forcefield) assigned by pdb2pqr
- **3_5pH_BSA.pdb**  
  Structure of Bovine Serum Albumin at pH 3.5 (CHARMM forcefield) assigned by pdb2pqr  

### mdp files
- **em.mdp**
- **nvt1.mdp**
- **nvt2.mdp**
- **npt.mdp**
- **nvtf_310K.mdp/nvtf_500K.mdp**  
  > Baseline: em -> nvt1 -> nvt2 -> npt -> nvtf_310K 
  > High temperatures: em -> nvt1 -> nvt2 -> npt -> nvtf_500K
  > Acidic pH: em -> nvt1 -> nvt2 -> npt -> nvtf_310K 


### LAMMPS .in file
- **nptDeform.in**  
  input file for shear simulation using LAMMPS
  > Shear: em -> nvt1 -> nvt2 -> npt -> convert to LAMMPS -> ptDeform.in

### Post analysis script
- **SS_func.py**  
  Analyze secondary structure from VMD timeline output (STRIDE)
