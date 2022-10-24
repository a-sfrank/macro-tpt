Copyright (c) 2022 Anna-Simone Frank and Susanna Röblitz with contributions from Markus Weber (in PCCA+)

If you use this code or parts of it, cite the following reference:

* Frank, A. S., Larripa, K., Ryu, H. and Röblitz S. (2022).
Macrophage phenotype transitions in a stochastic gene-regulatory network model. bioRxiv preprint:10.1101/2022.10.21.513139. Preprint available at https://biorxiv.org/cgi/content/short/2022.10.21.513139v1

* Frank, A. S., Sikorski, A., & Röblitz, S. (2022). 
Spectral clustering of Markov chain transition matrices with complex eigenvalues. arXiv preprint arXiv:2206.14537. Preprint available at https://doi.org/10.48550/arXiv.2206.14537.

* Röblitz, S., & Weber, M. (2013).
Fuzzy spectral clustering by PCCA+: application to Markov state models and data classification. Advances in Data Analysis and Classification, 7(2), 147-179.https://doi.org/10.1007/s11634-013-0134-6

MACRO-TPT is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

------------------------------------------------------------------------------------------
For questions or support contact

Anna-Simone Frank (asfrank88@gmail.com; anna-simone.frank@uib.no) or Susanna Röblitz (susanna.roblitz@uib.no)

------------------------------------------------------------------------------------------

Download the code at: https://github.com/a-sfrank/macro-tpt.git

------------------------------------------------------------------------------------------
# General information:
The code implements the PCCA+ method [1, 2] and performs Transition path theory analysis [3--5] 
on the coarse-grained Markov State Model (MSM) to determine transition statistics, path flows and transition probabilities.
The code is applied to a stochastic model representation (Chemical Master Equation) of a macrophage polarization model in [6] with the aim to identify multi-stable macropahge phenotypes and the transitions between them.

Code generates the results presented in article:
Frank, A. S., Larripa, K., Ryu, H. and Röblitz S. (2022).
*Macrophage phenotype transitions in a stochastic gene-regulatory network model*. bioRxiv preprint:10.1101/2022.10.21.513139. Preprint available at https://biorxiv.org/cgi/content/short/2022.10.21.513139v1

The code implementation is based on the following articles: 
1. Frank, A. S., Sikorski, A., & Röblitz, S. (2022). Spectral clustering of Markov chain transition matrices with complex eigenvalues. arXiv preprint arXiv:2206.14537. Preprint available at https://doi.org/10.48550/arXiv.2206.14537
2. Röblitz, S., & Weber, M. (2013). Fuzzy spectral clustering by PCCA+: application to Markov state models and data classification. Advances in Data Analysis and Classification, 7(2), 147-179. https://doi.org/10.1007/s11634-013-0134-6
3. Metzner, P., Schütte, C., & Vanden-Eijnden, E. (2009). Transition path theory for Markov jump processes. Multiscale Modeling & Simulation, 7(3), 1192-1219. https://doi.org/10.1137/070699500
4. Helfmann, L., Ribera Borrell, E., Schütte, C., & Koltai, P. (2020). Extending transition path theory: Periodically driven and finite-time dynamics. Journal of nonlinear science, 30(6), 3321-3366.https://doi.org/10.1007/s00332-020-09652-7
5. Noé, F., Schütte, C., Vanden-Eijnden, E., Reich, L., & Weikl, T. R. (2009). Constructing the equilibrium ensemble of folding pathways from short off-equilibrium simulations. Proceedings of the National Academy of Sciences, 106(45), 19011-19016. https://doi.org/10.1073/pnas.0905466106
6. Frank, A. S., Larripa, K., Ryu, H., Snodgrass, R. G., & Röblitz, S. (2021). Bifurcation and sensitivity analysis reveal key drivers of multistability in a model of macrophage polarization. Journal of Theoretical Biology, 509, 110511. https://doi.org/10.1016/j.jtbi.2020.110511
7. Chu, B. K., Tse, M. J., Sato, R. R., & Read, E. L. (2017). Markov State Models of gene regulatory networks. BMC systems biology, 11(14), 1-17. https://doi.org/10.1186/s12918-017-0394-4

# Included code files and their description:
|File names:		|	Description:
|---|---|
|Input.m:| 		 	File includes and specifies input parameters important for PCCA+ and TPT analysis|
|main_pcca_macrophages.m:| 	Main file that loads the stochastic macrophage model, pcca+ and TPT analysis|
|parameterCase.m:	| 	File contains the different model parameter cases: 1.	**Bistability:** Case 1-4; 2. **Tristability:** Case 5 with subcases 1-6|
|assemble Q.m:	|		File assembles the stochastic rate matrix from the macrophage polarization model described in Frank et al. (CME generation)|
|pcca.m:|			 	Main file for the pcca+ analysis. It depends on the following files (compute_subspace.m; fillA; index_search.m; main_nlscon.m; nlscon.m; objective.m; orthogon.m; opt_soft.m; problem_pcca-nlscon.m). For details on the pcca+ method and the code, see Röblitz et al.[2] and Frank et al.[1] as well as https://github.com/sroeblitz/cPCCA.git  (Main PCCA+ analysis file) |
|pikFullfnc.m:|			Function calculates the density peaks from the clustering result from the PCCA+ analysis.|
|findSet_2.m:| 			Function based on pikFullfnc.m and is used to define the macrophage phenotypes sets on the state space. These phenotype sets are needed to calculate and perform TPT. They are loaded into the TPT_cases.m file.|
|committor.m:|			Code implements committor function as described in article by Metzner et al.[3] (Part of TPT analysis)|
|prob2.m; prob3.m, prob4.m:|   Code implements Transition path theory (TPT): Calculates transition flows between sets and probabilities and	stationary probability of the sets (i.e., phenotypes) for bistability (index 2), for tristability, subcases 1-2	(index 3) and for tristability, subcases 3-6 (index 4) |
|probBi.m, probTri.m:	|	Calculates additional statistics for the TPT analysis: for bistability (index: Bi);  for tristability (index: Tri)|
|TPTCases.m:|			File calls for each parameter case the files and codes needed to perform the TPT analysis. (Main TPT analysis file)|
|coords.m |                      The function calculated the 2D coordinates  of an index node and is needed to visualize the transition flows.  |
|plot_Feff_bistability.m| Function plots the flow directions over the state space N^2 for the *bistable* cases. One needs to specify which of the bistable cases shall be plotted. This function calls the coords.m function. |
|plot_Feff.m| Function plots the flow directions over the state space N^2 for the *tristable* cases. One needs to specify the corresponding case (i.e., Case 5), and subcase. This function calls the coords.m function. |
        
# Setting Input parameters:
|File names:|			Description:|
|---|---|
|Input.m:| 			Change and adapt input parameters for pcca+ and TPT analsyis in Input.m file; Define the filepath to folder where the results of analysis shall be saved.	|
|TPTCases.m:| 	**For Bistability:** To vary flow/transition direction adapt start and end set in the file. For Case 1 (lines 37/38), for case 2 (lines  72/73), for case 3 (lines 109/110), for case 4 (lines 146/147). *Note:* Start set is A; End set is B; Set C is the complement to the sets A and B in the state space. Set A or B represent one of the four macrophage phenotypes (LL, LH, HL or HH). **For Tristability:** NO ADAPTATIONS NEEDED. The variation in flow directions is accounted for in the various subcases.|

# Code outputs:
#### 1. From PCCA+ analysis
|Output type:	|		Description:|
|---|---|
|Figures|				Fig 1: Cluster vs Crispness plot;	Fig 2-3(4): Chi Membership functions for bi-(tri-)stability; Fig 11-12(13): Partial stationary densities of cluster sets (i.e.,macrophage phenotype sets)|
|Data |				Eigenvalues, Crispness of cluster; statistical weights of the cluster, transition rate matrix Qc and coarse-grained transition matrix Pc|

#### 2. From TPT analysis
|Output type:|			Description:|
|---|---|
|Figures|				Fig 21/22: Forward and Backward committors;	Fig 31/32: Respectively, probability and normalized probability distribution of reactive trajectories (mR and mAB);	Fig 33:	   Trantition path flow directions between phenotype sets|
|Data | Test result of markov chain reversibility; TPT statistics (Transition rates and transition times (tAB,kAB1,kAB2,ZAB);	Stationary probability of states (or cluster) (piA, piB,piB1,piB2); Transition probabilities within  and between states (e.g., TAA and TAB)--not used in above article; Total transition flux from A and into B (Flux out of A, Flux into B);	Flow decompositions for coarse grained flux (e.g., Flux A-->B, or C-->B). *Note:* The Flux information is used to calculate the transition flow probabilities (done manually); For the interpretation of this information, see above articles. All of the above data and Figures are saved automatically in the corresponding folder, which is specified by the given filepath (see Input.m file).|

# To run the code:
Perform the following steps:
1. Define input parameters in Input.m and the flow directions for Bistability analysis in TPTCases.m;
2. Run the main file main_pcca_macrophages.m: You are prompted to
	- select parameter case you want to run (Choose between 1-5 cases)
	- select the numbers of clusters, based on crispness value (see Fig 1) (select 2 for bistability, or 3 for tristability cluster)
	- in the tristable case (Case 5) select the subcase you want to run (choose between subcases 1 and 6)

# Postprocessing of code output (i.e., data):
- The flux information is used to calculate transition flow probabilities manually.
- Similarily, the calculated transition probabilities could be used to calculate/approximate the coarse-grained transition matrix Pc (see appendix 
 in [7]); This is NOT NECESSARY for replicating the results, as the PCCA+ method gives out the exact matrix Pc.
