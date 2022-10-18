README

Implements analysis in article:
"Stochastic framework reveals phenotype transitions in a model of macrophage polarization" by Frank Anna-Simone, Kamila Larripa, Hwayeon Ryu 
and Susanna Röblitz
==============================================================================================================================================
General information:
==============================================================================================================================================
The code implements the cPCCA+ method (Frank et. al; Röeblitz et al;) and performs Transition path theory analysis (Metzner et al.) 
on the coarse-grained Markov State Model (MSM) to determine transition statistics, path flows and transition probabilities.
The code is applied to a stochastic model representation (Chemical Master Equation) of a macrophage polarization model (by Frank et al.)
 with the aim to identify multi-stable macropahge phenotypes and the transitions between them.

It is based on the following articles: (Add articles here later)
- Frank et al
- Röblitz et al.
- Metner paper
- cPCCA paper
- Chu et al.
- etc.

Contributer: Anna-Simone Frank, Susanna Röblitz

Maintainer:

==============================================================================================================================================
Included code files  and their description:
==============================================================================================================================================
File names:			Description:
----------------------------------------------------------------------------------------------------------------------------------------------
Input.m: 		 	File includes input parameters important for PCCA+ and TPT analysis
main_pcca_macrophages.m: 	Main file that loads the stochastic macrophage model, pcca+ and TPT analysis
parameterCase.m:	 	File contains the different model parameter cases:
			 	Cases of bistability: Case 1-4
			 	cases of tristability: Case 5 with subcases 1-6
assemble Q.m:			File assembles the stochastic rate matrix from the macrophage polarization model described in Frank et al. (CME)
pcca.m:			 	Main file for the pcca+ analysis. It depends on the following files (compute_subspace.m; fillA; index_search.m;
			 	main_nlscon.m; nlscon.m; objective.m; orthogon.m; opt_soft.m; problem_pcca-nlscon.m). For details on the pcca+ 
				and the code, see Röblitz et al and Frank et al. (main PCCA+ analysis file)

pikFullfnc.m:			Function calculates the density peaks from the clustering result from the PCCA+ analysis.
findSet_2.m: 			Function based on pikFullfnc.m and is used to define the macrophage phenotypes sets on the state space.
committor.m:			Code implements committor function as described in article by Metzner et al. (part of TPT analysis)

prob2.m; prob3.m, prob4.m:      Code implements Transition path theory (TPT): Calculates transition flows between sets and probabilities and
				stationary probability of the sets (i.e., phenotypes) for bistability (index 2), for tristability, subcases 1-2
				(index 3) and for tristability, subcases 3-6 (index 4) 
probBi.m, probTri.m:		Calculates additional statistics for the TPT analysis and plots the flow directions over the state space N^2:
				for bistability (index: Bi), for tristability (index: Tri)
TPTCases.m:			File calls for each parameter case the files and codes needed to perform the TPT analysis. 
				(Main TPT analysis file)

==============================================================================================================================================
Input parameters:
==============================================================================================================================================	
File names:			Description:
----------------------------------------------------------------------------------------------------------------------------------------------
Input.m: 			Change and adapt input parameters for pcca+ and TPT analsyis in Input.m file; Define the filepath to folder 
				where the results of analysis shall be saved.	

TPTCases.m: 			For Bistability: 
				- To vary flow/transition direction vary start and end set in the file. For Case 1 (lines 37/38), 
				for case 2 (lines  72/73), for case 3 (lines 109/110), for case 4 (lines 146/147).
				Note:
				- Start set is A; End set is B; C is everything that is the complement of the sets A and B to the rest of the 
				state space.
				- Dependend on the flow direction set A or B are macrophage phenotypes (LL, LH, HL or HH)

				For Tristability:
				- NO ADAPTIONS NEEDED. The variation in flow directions is accounted for in the various subcases.

==============================================================================================================================================
Code outputs:
==============================================================================================================================================	
1. From PCCA+ analysis
Output type:			Description:
----------------------------------------------------------------------------------------------------------------------------------------------
Figures				Fig 1: Cluster vs Crispness plot
				Fig 2-3(4): Chi Membership functions for bi-(tri-)stability
				Fig 11-12(13): Partial stationary densities of cluster sets (i.e.,macrophage phenotype sets)

Data 				Eigenvalues, Crispness of cluster; statistical weigthes of the cluster and coarse-grained transition matrix Pc

2. From TPT analysis
Output type:			Description:
----------------------------------------------------------------------------------------------------------------------------------------------
Figures				Fig 21/22: Forward and Backward committors
				Fig 31/32: Respectively, probability and normalized probability distribution of reactive trajectories (mR and mAB)
				Fig 33:	   Trantition path flow directions between phenotype sets

Data  				Test result of markov chain reversibility; TPT statistics (Transition rates and transition times (tAB,kAB1,kAB2,ZAB),
				Stationary probability of states (or cluster) (piA, piB,piB1,piB2), Transition probabilities within  and between  
				states (e.g., TAA and TAB)--not used; Total transition flux from A and into B (Flux out of A, Flux into B),
				Flow decompositions for coarse grained flux (e.g., Flux A-->B, or C-->B)
				Note: The Flux information is used to calculate the transition flow probabilities (done manually); For the 
				interpretability of this information, see above articles.

All of the above data and Figures are saved automatically in the corresponding folder, which is specified by the given filepath (see Input.m file).

==============================================================================================================================================
How to run the code:
==============================================================================================================================================
Perform the following steps:
1. Define input parameters in Input.m and the flow directions for Bistability analysis in TPTCases.m;
2. Run the main file main_pcca_macrophages.m: You are prompted to
	- select parameter case you want to run (Choose between 1-5 cases)
	- select the numbers of clusters, based on crispness value (see Fig 1) (select 2 for bistability, or 3 for tristability cluster)
	- in the tristable case (Case 5) select the subcase you want to run (choose between subcases 1 and 6)

Postprocessing of data:
- The flux information is used to calculate transition flow probabilities manually.
- Similarily, the calculated transition probabilities could be used to calculate/approximate the coarse-grained transition matrix Pc (see appendix 
 in Chu et al.); For this article this is NOT DONE, as the PCCA gives out the exact matrix PC;

==============================================================================================================================================
OPEN issues: (need to be addressed in the AIM^2 Workshop week)
==============================================================================================================================================
There are still open issues that have not yet been fixed or determined: (the list below might not be complete)

1. Best size of state-space N needs to be determined. (Input.m file)
2. Results depend on the value of the parameter PercTh (see Input.m file) and so does the state-space N size (I think)  
3. Flow direcitons in the cases 1-4 (Bistability) have to be determined (in TPTCases.m file), they are not yet fixed for the bistable cases.
4  I deemed in Case 5--subcase 2 as uninteresting. Everybody shall check if they agree with this choice.
