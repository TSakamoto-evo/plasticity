###################################################
Codes for simulations in
Sakamoto and Innan (2023)

We updated the model during the revision and the current model is different from the old model which is used in the preprint version uploaded in June 2022.
C++ scripts of the old model are in "old" folder.
###################################################

Contents
1. Folder "simulation_main"
C++ scripts for simulations of Wright-Fisher model that was used to generate results in Table 1 and Table 2 are provided.
Free recombination is assumed.

(How-to-use)
	(i) Compile by the following command.
		"g++ *.cpp -std=c++11 -O3 -o XXX.out"
	(ii) Run the executive file by the following command. Each xi denotes command line arguments.
		"./XXX.out x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19"
	(iii) Instead, for the parameter set in Table 1, we can run a simulation by the following command:
		"python3 run_each_para.py X" where X is the index for a parameter set. See the file "parameter_list_raw.csv" for details.

(Command line arguments)
	x1: int type. An arbitrary index that is used for distinguishing simulation runs.
	x2: int type. Haploid population size, N.
	x3: int type. The number of developmental steps, tau_v.
	x4: int type. The number of generations between switches of environments, I_T.
	x5: double type. Degeneration rate of gene expression, alpha.
	x6: double type. Maximum mutation impact for each g_i, gamma_g.
	x7: double type. Maximum mutation impact for each b_ij, gamma_b.
	x8: double type. Maximum mutation impact for each c_i, gamma_c.
	x9: double type. Strength of selection, sigma.
	x10: double type. Criteria for the evolution of plasticity. In all cases, 0.95 is used.
	x11-18: int type. Optimal expression patterns in two environments. Only following four values are allowed: 0, 1, 2, -2.
		0 represents that optimal expressions are 0 in both environments.
		1 represents that optimal expressions are 1/alpha in both environments.
		2 represents that optimal expressions are 1/alpha in E1 and 0 in E2.
		-2 represents that optimal expressions are 0 in E1 and 1/alpha in E2.
		For example, (x11, x12, x13, x14, x15, x16, x17, x18) = (0, 0, 1, 1, 2, 2, -2, -2) when alpha = 0.2 means 
		X1 = (0, 0, 5, 5, 5, 5, 0, 0) and X2 = (0, 0, 5, 5, 0, 0, 5, 5).
	x19: int type. An arbitrary index that is used for distinguishing parameter sets.

(Expected output)
Time until plasticity is evolved, T_p, is written in file "regi_time.txt" in append mode in the following format.
	"simulation_index(x1)	parameter_index(x19)	whether_evolution_is_successful(1)_or_not(0)	time_until_plasticity_is_evolved(T_p, scaled by N)"

The results of development in hypothetical intermediate environments "regi_env_grad.txt" in append mode in a following format.
	"simulation_index(x1)	parameter_index(x19)	whether_evolution_is_successful(1)_or_not(0)	# of environments in which w1 > criteria	# of environments in which w2 > criteria	# of environments in which both w1 and w2 are below the criteria
 

2. Folder "simulation_linked"
C++ scripts for simulations of Wright-Fisher model that was used to generate Table S1 are provided.
Only the evolution of cis-regulatory is considered, and all cis-regulatory regions of a certain gene are assumed to be completely linked.
Usage is the same as the folder "simulation_main".

3. Folder "simulation_visualization"
C++ scripts and R scripts that were used to generate figures 1, 2, S1, and S2 are provided.
Core of the C++ simulation is the same as "simulation_main", but it generates more detailed outputs.

(c++ output files)
	"regi_genotypes.txt": Return the genotype of individuals.
	"answers.txt": Return parameter values.
	"fitness.txt": Return mean fitnesses (w1 and w2) in the two environments.

R scripts are used to visualize C++ output.
	"plot_fitness.R": plot mean fitnesses (w1 and w2).
	"plot_genotypes.R": plot genotype and phenotype. Specify the time when individuals are sampled in the file.
	"plot_development.R": plot genotype and phenotype. Specify the time when individuals are sampled in the file.
	

(How-to-use)
	(i) Run C++ simulations in the same way as the above cases.
	(ii) Run R scripts.

