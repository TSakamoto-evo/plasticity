#ifndef GENOTYPE
#define GENOTYPE

#include <vector>
#include <cmath>
#include <iostream>

class Genotype{
private:
  std::vector<double> genotype;
  std::vector<std::vector<double>> interactions;
  std::vector<double> env_effects;
  int num_g;

  int ind_num;
  std::vector<std::vector<double>> phenotypes;
  int adult_step;

public:
  Genotype(const int input_ind_num, const std::vector<double>& input_genotype,
    const std::vector<std::vector<double>>& input_interactions,
    const std::vector<double>& input_env_effects);
  void set_genotype(const int input_ind_num,
    const std::vector<double>& input_genotype,
    const std::vector<std::vector<double>>& input_interactions,
    const std::vector<double>& input_env_effects);
  void development(const double alpha, const int round, const double env_cue,
    const int sensor_period);
  double ret_fitness(const double alpha, const int round, const double env_cue,
    const int sensor_period, const std::vector<double>& s_vec, const double deno);
  void add_mut(const int index, const double delta, const double max_g);
  void return_genotype(std::vector<double>& ret) const{ ret = genotype; };
  void return_interactions(std::vector<std::vector<double>>& ret) const
    { ret = interactions; };
  void return_env_effects(std::vector<double>& ret) const{ ret = env_effects; };
  void return_final(const std::vector<double>& ini_g,
    std::vector<std::vector<double>>& ret, const double alpha, const int round) const;
  void return_final_grad(const double env_cue,
    std::vector<std::vector<double>>& ret, const double alpha, const int round,
    const int sensor_period) const;
  bool check_identical(const std::vector<double>& input_genotype,
    const std::vector<std::vector<double>>& input_interactions,
    const std::vector<double>& input_env_effects) const;
  void add_ind(){ ind_num++; };
  void reduce_ind(){ ind_num--; };
  void set_ind_num(const int num){ ind_num = num; };
  int ret_ind_num() const{ return(ind_num); };
  bool fluct(const int wait, const std::vector<double>& env_cues,
    const double alpha, const int round, const int sensor_period) const;
};



#endif
