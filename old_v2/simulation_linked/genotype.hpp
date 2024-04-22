#ifndef GENOTYPE
#define GENOTYPE

#include <vector>
#include <cmath>
#include <iostream>

class Genotype{
private:
  int ind_num;
  std::vector<double> genotype;
  std::vector<std::vector<double>> interactions;
  std::vector<double> env_effects;

  std::vector<int> mut_allele;

  std::vector<double> phenotype;
  int num_g;

  double fitness;

public:
  Genotype(const std::vector<double> input_genotype,
    const std::vector<std::vector<double>> input_interactions,
    const std::vector<double> input_env_effects,
    const std::vector<int> input_mut_allele,
    const int input_ind_num);
  void development(const double alpha, const int round, const double env_cue,
    const int sensor_period);
  void set_fitness(const std::vector<double> s_vec, const double sig);
  double return_geno_value(const int index, int& allele_index) const;
  void set_mutation(const std::vector<int>& locus_index,
    const std::vector<int>& index_list,
    const std::vector<std::vector<double>>& effect_list,
    const int input_ind_num);
  int return_size() const{ return(ind_num); };
  void reduce_ind(){ ind_num--; };
  std::vector<double> add_mut(const int index, const double delta, const int allele_index,
    const double max_g);
  double return_fit_size() const{ return(ind_num * fitness); };
  double return_fitness() const{ return(fitness); };
  void return_phenotype(std::vector<double>& ret) const{ ret = phenotype; };
  void return_genotype(std::vector<double>& ret) const{ ret = genotype; };
  void return_env_effects(std::vector<double>& ret) const{ ret = env_effects; };
  void return_interactions(std::vector<std::vector<double>>& ret) const;
  void return_final(const std::vector<double>& ini_g,
    std::vector<double>& ret, const double alpha, const int round);
  void return_final_grad(const double env_cue,
    std::vector<double>& ret, const double alpha, const int round,
    const int sensor_period);
};

#endif
