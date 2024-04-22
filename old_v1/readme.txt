###################################################
Codes for simulations in
Sakamoto and Innan (2022)
###################################################

Requirements: C++11 (tested with GCC 7.3.0 and clang 13.1.6), R (tested with version 4.0.5)
No non-standard hardware is required.
Installation time: If R and C++ is already installed, no appreciable time is needed.
Tested environments: MacOS (version 12.3), CentOS (version 7)

Contents
1. Folder "simulation_main"
C++ scripts for simulations of Wright-Fisher model that was used to generate results in Table 1 and Table S1 are provided.
Free recombination is assumed.

(How-to-use)
	(i) Compile by following command.
		"g++ *.cpp -std=c++11 -O3 -o XXX.out"
	(ii) Run the executive file by a following command. Each xi denotes command line arguments.
		"./XXX.out x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19"

(Command line arguments)
	x1: int type. Arbitrary index that is used for distinguishing simulation runs.
	x2: int type. Haploid population size, N.
	x3: int type. The number of developmental steps, v.
	x4: int type. The number of generations between switches of environments, T.
	x5: double type. Degeneration rate of gene expression, alpha.
	x6: double type. Maximum mutation effect of each g_i, mu_g.
	x7: double type. Maximum mutation effect of each b_ij, mu_b.
	x8: double type. Maximum mutation effect of each c_i, mu_c.
	x9: double type. Strength of selection, sigma.
	x10: double type. Criteria for evolution of plasticity. In all cases, 0.95 is used.
	x11-18: int type. Optimal expression patterns in two environments. Only following four values are allowed: 0, 1, 2, -2.
		0 represents that optimal expressions are 0 in both environments.
		1 represents that optimal expressions are 1/alpha in both environments.
		2 represents that optimal expressions are 1/alpha in E1 and 0 in E2.
		-2 represents that optimal expressions are 0 in E1 and 1/alpha in E2.
		For example, (x11, x12, x13, x14, x15, x16, x17, x18) = (0, 0, 1, 1, 2, 2, -2, -2) when alpha = 0.2 means 
		X1 = (0, 0, 5, 5, 5, 5, 0, 0) and X2 = (0, 0, 5, 5, 0, 0, 5, 5).
	x19: int type. Arbitrary index that is used for distinguishing parameter sets.

(Expected output)
Time until plasticity is evolved, T_p, is written in file "regi_time.txt" in append mode in a following format.
	"simulation_index(x1)	parameter_index(x19)	whether_evolution_is_successful(1)_or_not(0)	time_until_plasticity_is_evolved(T_p)"

(Run time)
Typical maximum running time of one run is ~8 hours in a normal computer within investigated parameter space.
 

2. Folder "simulation_linkage2"
C++ scripts for simulations of Wright-Fisher model that was used to generate Table S2 are provided.
Only the evolution of cis-regulatory is considered, and all cis-regulatory regions on a certain gene are assumed to be completely linked.
Usage is the same as the folder "simulation_main".

3. Folder "simulation_visualization"
C++ scripts and R scripts that were used to generate figures 1, 2, and S1 are provided.
Core of the C++ simulation is the same as "simulation_main", but it generates more detailed output.

(c++ output files)
	"regi_genotypes.txt": Return the genotype of a randomly-sampled individual.
	"answers.txt": Return optimal expression levels in two environments.
	"fitness.txt": Return mean fitnesses (w1 and w2) in the two environments.

R scripts are used to visualize C++ output.
Explanation for each R script is as follows:
	"plot_fitness.R": plot mean fitnesses (w1 and w2).
	"plot_genotypes_one.R": plot genotype and phenotype. Specify which row of "regi_genotype.txt" is used in the file.
	"plot_dev_process_one_time.R": plot genotype and phenotype. Specify which row of "regi_genotype.txt" is used in the file.
	

(How-to-use)
	(i) Set appropriate parameters in main.cpp file.
	(ii) Compile C++ code by following command.
		"g++ *.cpp -std=c++11 -O3 -o XXX.out"
	(iii) Run the executive file by a following command.
		"./XXX.out"
	(iv) Run R scripts.

(Run time)
Typical maximum running time of one run is ~8 hours in a normal computer within investigated parameter space.

