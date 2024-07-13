
Readme for TiRNA package by Tan-group at Wuhan University

TiRNA: a coarse-grained method with temperature and ion effects for RNA structure folding and prediction 

Overview of the TiRNA
The model has three main features:
(1) A three-bead coarse-grained representation: P, C4' and N1/N9;
(2) A coarse-grained force field involving effects of temperature and ions; 
(3) Replica-exchange Monte Carlo/Monte Carlo simulated annealing for structure sampling.

The package of the model contains two modules:  
(1) RNA 3D structure prediction: to predict 3D structures of RNAs from sequences and certain monovalent/divalent ion conditions or from secondary structures; 
(2) RNA thermal stability prediction: to predict melting temperatures and thermal unfolding pathway from sequences and certain monovalent/divalent ion conditions.  

Machine requirements:
(1) Linux operating systems;
(2) >=8 CPU cores with 2 threads;
(3) GCC version (>= 7.5) supporting c++11 standard;
(4) Python version (>= 3.11.5).


Readme for running TiRNA: 
__________________________________________________________________

1, RNA 3D structure and thermal stability prediction from sequences.
(Run in folder1)
    
(1) Inputting an RNA sequence in seq.dat

(2) Inputting simulation configuration parameters in the config.dat
*configuration parameters*
Sampling_type 1
Folding_steps 750000   	   
Optimizing_steps 500000     
C_Na 1000   	           
C_Mg 0			    
N_cout 10 
                   
Notes:

Sampling_type: 1 and 0 corresponds to the sampling by replica-exchange Monte Carlo and the sampling by Monte Carlo simulated annealing; 
Folding_steps: Number of steps for replica-exchange Monte Carlo or Monte Carlo annealing simulations in structure folding progress;
Optimizing_steps: Number of steps in structure optimizing progress;
C_Na: Concentration of Na (in mM);
C_Mg: Concentration of Mg (in mM);
N_cout: Number of predicted 3D structures.
The above shown parameters are our suggested ones.

(3) bash run.sh
The predicted structures are given in a new directory "results", see "example1" for an example.
(a) The results are placed in the 'folder1'; 
(b) Folding_trajectory contains folding trajectories at different temperatures;
    CG_structures contains the predicted coarse-grained top 3D structures;
    Secondary_structure contains the predicted top-N secondary structures in dot-bracket form;
    All_atom_structure contains All_atom_structure corresponding top-N CG structures;
    Thermal_Stability contains the fractions of different states (Folded, Unfolded, Intermediate) at different temperatures recorded in thermal_stabiliy.dat' and the secondary structures in dot-bracket form at different temperatures.

Notes: 	
(a) If you want to predict RNA 3D structures more accurately, the suggested step numbers should be >750000 (500000 is the minimum) for the folding steps and >500000 (100000 is the minimum) for the optimizing steps;
(b) If you want to predict RNA thermal stability more accurately, the suggested step number should be >4000000 (2000000 is the minimum) for the folding steps; 
(c)  If you want to obtain the melting temperatures, please fit fractions of Folded or Unfolded states to a two-state model according to the references [1-3]. 



2, RNA 3D structure and thermal stability prediction with given the secondary structures.
(Run in folder2)

(1) Inputting an RNA sequence and its secondary structure (dot-bracket form) in seq.dat of 'folder2';
(2) The other operation processes are the same as (A) for inputting only sequence.


3, Further refinement is required for the rebuilt all-atom structures.

To remove possible steric clashes and chain breaks of the rebuilt all-atom structures by TiRNA, a structure refinement can be performed through the method of QRNAS (https://github.com/sunandan-mukherjee/QRNAS.git) [4].

__________________________________________________________________

If you have any questions about TiRNA, please contact us by the email: zjtan@whu.edu.cn.

References:
[1] Wang X, Lou E, Yu S, Tan YL, Shi YZ, & Tan ZJ. 2024. TiRNA: a coarse-grained method with temperature and ion effects for RNA structure folding and prediction. In preparation.
[2] Wang X, Tan YL, Yu S, Shi YZ, & Tan ZJ. 2023. Predicting 3D structures and stabilities for complex RNA pseudoknots in ion solutions. Biophys J. 122, 1503-1516.
[3] Shi YZ, Wang FH, Wu YY, & Tan ZJ. 2014. A coarse-grained model with implicit salt for RNAs: Predicting 3D structure, stability and salt effect. J Chem Phys. 2014, 141, 105102.
[4] Stasiewicz J, Mukherjee S, Nithin C, & Bujnicki JM. 2019. QRNAS: Software tool for refinement of nucleic acid structures. BMC Struct Biol. 19, 5.
