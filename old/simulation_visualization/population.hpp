#ifndef POPULATION
#define POPULATION

#include "genotype.hpp"
#include "mutation.hpp"
#include "parameters.hpp"
#include <vector>
#include <random>
#include <cmath>
#include <iostream>

class Population{
private:
  std::vector<Genotype> pop;
  Mutation mutation;
  Parameters para;
  std::mt19937 mt;
  int gen;

  std::poisson_distribution<> mut_num;
  std::discrete_distribution<> choose_type;
  std::uniform_int_distribution<> choose_locus;
  std::uniform_real_distribution<> d_mu_g;
  std::uniform_real_distribution<> d_mu_b;
  std::uniform_real_distribution<> d_mu_c;

public:
  Population(const Parameters input_para);
  void development(const double env_cue, const std::vector<double> s_vec);
  void selection();
  int return_total_geno() const;
  void gen_mutation();
  void one_generation(const double env_cue, const std::vector<double> s_vec);
  void selection_many_segregation();
  std::vector<double> calculate_mean_fitness(const std::vector<double>& env_cues,
    const std::vector<std::vector<double>> s_vecs);
  void return_ind(std::vector<double>& ret_genotypes,
    std::vector<std::vector<double>>& ret_interactions,
    std::vector<double>& ret_env_effects,
    std::vector<double>& ret_adult1, std::vector<double>& ret_adult2,
    const double s1_env, const double s2_env);
  void return_dev_process(const double s1_env, const double s2_env,
    std::vector<std::vector<double>>& ret_dev1,
    std::vector<std::vector<double>>& ret_dev2);
  std::vector<std::vector<double>> random_genotype(const int trial,
    const double s1_env, const double s2_env);
  std::vector<std::vector<double>> grad_env(const std::vector<double>& grad_cues);
};

#endif
