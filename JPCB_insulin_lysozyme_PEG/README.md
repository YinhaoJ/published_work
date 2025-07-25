# PEGylation of Insulin and Lysozyme To Stabilize against Thermal Denaturation: A Molecular Dynamics Simulation Study
This repository contains input files, scripts, and protein structures used in the simulations for the following publication:

*PEGylation of Insulin and Lysozyme To Stabilize against Thermal Denaturation: A Molecular Dynamics Simulation Study*


---

## Content

### Protein Structure
- **Lysozyme_K33_PEG75_nvtq.pdb**  
  Equilibrated Lysozyme_K33_PEG50 structure 
- **Insulin_PEG50_nvtq.pdb**  
  Equilibrated Insulin_BK29_PEG50 structure  

### mdp files
- **em.mdp**
- **nvt1.mdp**
- **nvt2.mdp**
- **npt.mdp**
- **nvtf.mdp**  
- **nvtq.mdp**
> for non-PEGylated em -> nvt1 -> nvt2 -> npt -> nvtf (300K or 500K) 
> for PEGylated em -> nvt1 -> nvt2 -> npt -> nvtq -> nvtf (300K or 500K)

### Post analysis script
- **Script.ipynb**  
  **Get native contact from crystal structure PDB**
  **Convert secondary structure data (VMD timeline) to dataframe and plot**

