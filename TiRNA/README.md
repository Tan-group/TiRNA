#TiRNA: a coarse-grained method with temperature and ion effects for RNA structure folding and prediction

##Overview
================================================================

TiRNA is a coarse-grained computational method developed by Tan-group at Wuhan University for RNA structure modeling. The model incorporates the effects of temperature and ions to provide predictions for 3D structures and thermal stability of RNAs in monovalent/divalent ion solutions.

###Key Features:
(1) Three-bead coarse-grained representation (P, C4', and N1/N9);
(2) Force field accounting for temperature and ion effects;
(3) Replica-exchange Monte Carlo & Monte Carlo simulated annealing for structure sampling;

###Package Modules:
(1) RNA 3D Structure Prediction: Predicting 3D structures from sequences or secondary structures at given ion conditions;
(2) RNA Thermal Stability Prediction: Predicting structures at various temperatures which can be used for analyzing melting temperatures and thermally unfolding pathways.



##System Requirements
================================================================
- **Operating System**: Linux;
- **CPU**: ≥8 cores with 2 threads;
- **Compiler**: GCC ≥7.5 (supporting C++11 standard);
- **Python**: Version ≥3.11.5;
- **Python Packages**: Biopython ≥1.80 (pip install biopython).



##Installation
================================================================

1. Download the TiRNA package
git clone https://github.com/Tan-group/TiRNA.git


2. Verify dependencies
gcc --version  # Should be ≥7.5;
python --version  # Should be ≥3.11.5;
python -c "import Bio; print(f'Biopython version: {Bio.__version__}')" #≥1.70 .


Usage Examples
================================================================

1. Prediction for RNA 3D Structure and Thermal Stability from Sequences 
Location: Run in from-sequence/

Steps:

(1) Prepare sequence file (seq.dat)  
# Input your RNA sequence in seq.dat
# Example:
GGCGAUGUCCAGCAGAUACACGUCGUUCGCACC

(2) Configure parameters (config.dat)
Sampling_type 1
Folding_steps 750000
Optimizing_steps 500000
C_Na 1000
C_Mg 0
N_cout 10

(3) Run the simulation with TiRNA
bash run.sh


2. RNA 3D Structure Prediction with Given Secondary Structure
Location: Run in from-2D-structure/

Steps:

(1) Prepare sequence and secondary structure (seq.dat)
# Format: sequence followed by secondary structure in dot-bracket notation
# Example:
GGAGGAAGGAGCCUCC
(((((......)))))
(2) Configure parameters (config.dat) - same as those in from-sequence
(3) Run the simulation with TiRNA
bash run.sh


3. RNA 3D Structure Prediction from Initial Structures

3A. From Initial Structure in PDB Format
Location: Run in from-3D-structure/PDB
Steps:
(1) Place your initial PDB structure file in the directory;
(2) Configure parameters in config.dat;
(3) Run with bash run.sh.


3B. From 3D Structure in CIF Format 
Location: Run in from-3D-structure/CIF
Steps:
(1) Place your initial CIF structure file in the directory;
(2) Configure parameters in config.dat;
(3) Run with bash run.sh.
    
    
    
    

Configuration Parameters
================================================================

## Configuration Parameters

| Parameter         | Description |
|-------------------|-------------|
| `Sampling_type`   | `1` for replica-exchange Monte Carlo (REMC),
                      `0` for Monte Carlo simulated annealing (MCSA) |
| `Folding_steps`   | Number of steps for structure folding |
| `Optimizing_steps`| Number of steps for structure optimization |
| `C_Na`            | Na⁺ concentration in mM |
| `C_Mg`            | Mg²⁺ concentration in mM |
| `N_cout`          | Number of output predicted 3D structures |

### Performance Recommendations:

#### For Replica-Exchange Monte Carlo (REMC, Sampling_type = 1):
- **Folding_steps**: ≥500,000 steps recommended;
- **Optimizing_steps**: ≥100,000 steps recommended.

#### For Monte Carlo Simulated Annealing (MCSA, Sampling_type = 0):
- **Folding_steps**: ≥2,000,000 steps recommended;
- **Optimizing_steps**: ≥100,000 steps recommended.

#### Application-specific recommendations
- For more accurate **3D structure prediction**: 
  - REMC: Use ≥750,000 folding steps and ≥500,000 optimizing steps;
  - MCSA: Use ≥3,000,000 folding steps and ≥500,000 optimizing steps.
- For more accurate **thermal stability prediction**: 
  - REMC: Use ≥4,000,000 folding steps;
  - MCSA: Use ≥20,000,000 folding steps.
    
    
    

Output Files
================================================================

After successful execution, results are saved in the results/ directory:
| Directory              | Content Description |
|------------------------|---------------------|
| `Folding_trajectory/` | Folding trajectories at different temperatures|
| `CG_structures/` | Predicted top-N 3D CG structures |
| `Secondary_structure/` | Predicted top-N secondary structures in dot-bracket notation |
| `All_atom_structure/`  | All-atom structures corresponding to top-N CG structures |
| `Thermal_Stability/`   | Data for thermal stability analysis |




Post-Processing
================================================================

To refine the predicted all-atom structures and remove potential steric clashes or chain breaks, use QRNAS:
# Install QRNAS from: https://github.com/sunandan-mukherjee/QRNAS.git




Support and Contact
================================================================

For questions or issues regarding TiRNA, please contact: zjtan@whu.edu.cn.



References
================================================================

[1] Wang X, Lou E, Yu S, Tan YL, Shi YZ, & Tan ZJ. 2025. TiRNA: a coarse-grained method with temperature and ion effects for RNA structure folding and prediction. In preparation.
[2] Wang X, Tan YL, Yu S, Shi YZ, & Tan ZJ. 2023. Predicting 3D structures and stabilities for complex RNA pseudoknots in ion solutions. Biophys J. 122, 1503-1516.
[3] Shi YZ, Wang FH, Wu YY, & Tan ZJ. 2014. A coarse-grained model with implicit salt for RNAs: Predicting 3D structure, stability and salt effect. J Chem Phys. 141, 105102.
[4] Stasiewicz J, Mukherjee S, Nithin C, & Bujnicki JM. 2019. QRNAS: Software tool for refinement of nucleic acid structures. BMC Struct Biol. 19, 5.

TiRNA Package - Tan Group, Wuhan University

