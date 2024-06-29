********** Readme for TiRNA software ******** By Tan-group at Wuhan University

# TiRNA
Predicting RNA 3D structure and thermal stability from the squeuence
Contributed by Xunxun Wang, Ya-Lan Tan, Ya-Zhou Shi and Zhi-Jie Tan 
This is free software for non-commercial users. If you have any questions, please contact Zhi-Jie Tan via email zjtan@whu.edu.cn


# Overview of the TiRNA
TiRNA is a coarse-grained model for simulating RNA structure folding and predicting with effects of temperatures and ios. 
The model has three main characteristics:
(1) A three beads coarse-grained representation: P, C4' and N1/N9.
(2) A coarse-grained force field involving temperature and ion effects. 
(3) Replica-exchange Monte Carlo or Monte Carlo annealing simulations for conformational sampling. 

The software contains two modules:  
(1) "RNA 3D Structure prediction": To predict 3D structures of RNAs only given the sequence and ion conditions or with given secondary structure constraints; 
(2) "RNA thermal stability prediction": To predict the melting temperature and theraml unfolding pathway only given the sequence and ion conditions.  


## Machine requirements:
Linux operating systems
cores >=8 cores and CPU should have each core paird with 2 threads
Gcc version (>= 7.5) supporting c++11 standard
Python version (>= 3.11.5)

### Usage:

***************************************************************************
(A) RNA 3D structure and thermal stability prediction from the sequence:

    (1) Inputing RNA sequence in seq.dat
    (2) Setting the simulation configuration parameters in the config.dat
    **** configuration parameters *****
    Sampling_type 1
    Folding_steps 750000   	   
    Optimizing_steps 10000     
    CNa 1000   	           
    CMg 0			    
    Ncout 10 
                   
    Explanation:
       Sampling_type 1: Sampling 1 means simulation by replica exchange Monte Carlo and Sampling 0 means simulations by Monte Carlo annealling 
       Folding_steps: Number of steps for replica exchange Monte Carlo or Monte Carlo annealing simulations in folding progress (suggestio)
       Optimizing_steps: Number of steps in optimizing progress
       CNa 10000: Concentration of Na (unit:mM)
       CMg 0: Concentration of Mg     (unit:mM)
       Ncout 10: Number of predicted 3D structures
    (3) bash run.sh
   
The predicted structures are given in a new directory "results", see "example1" for an example.
- (a) The results are placed in the 'folder1'; 
- (b) 'Folding_trajectory' contains folding trajectory;
      'CG_structures' contains the predicted coarse-grained structures;
      'Secondary_structure' contains the secondary structure (dot and bracket form) corresponding CG structures;
      'All_atom_structure' contains the All_atom_structure corresponding CG structures;
      'Thermal_Stability' contains the the fractions of different states (Folded, Unfolded, Intermediate) at different tempeartures, and recorded in 'thermal_stabiliy.dat' and the secondary structures (dot and bracket form) with different conformations at different temperatures.
*****Note: 	
	If you want to predict RNA 3D structure more accurately, the suggested steps for the folding_steps is >750000 (must be > 500000),optimizing_steps is >500000 (must be > 100000)
	If you want to predict RNA thermal stability more accurately, the suggested steps for the folding_steps is >4000000 (must be > 2000000)
	If you want to obtain the melting temperature, please refer to the reference [1,2]. 


************************************************************************************
(B) RNA 3D structure and thermal stability prediction with given the secondary structure constraints:
	
	The operation process is the same as (A).

## Further refinement is required for the rebuilt all-atom structures.
To remove possible steric clashes and chain breaks of the rebuilt all-atom structures,  a structure 
refinement  can be performed for the rebuilt all-atom structures by FebRNA through the method 
of QRNAS (https://github.com/sunandan-mukherjee/QRNAS.git) [3].
