#ifndef MUTATION
#define MUTATION

#include <vector>
#include <unordered_map>
#include <iostream>

class Mutation{
private:
  std::vector<std::unordered_map<int, double>> mut_genotype;
  std::vector<std::vector<std::unordered_map<int, double>>> mut_interactions;
  std::vector<std::unordered_map<int, double>> mut_env_effects;

  std::vector<int> index_genotype;
  std::vector<std::vector<int>> index_interactions;
  std::vector<int> index_env_effects;

  std::vector<int> num_genotype;
  std::vector<std::vector<int>> num_interactions;
  std::vector<int> num_env_effects;
  int num_g;

public:
  Mutation(const int input_num_g);
  void add_mutation(const int index, const double effect);
  double return_mut_effect(const int index, const int key) const;
  int allele_num(const int index) const;
  std::vector<int> return_allele_index(const int index,
    std::vector<double>& effects) const;
  void erase(const int locus_index, const int allele_index);
  int return_index(const int index) const;
  void check_status() const;

  void ret_mut_genotype(const int index,
    std::vector<int>& ret_key, std::vector<double>& ret_value) const;
  void ret_mut_interaction(const int index1, const int index2,
    std::vector<int>& ret_key, std::vector<double>& ret_value) const;
  void ret_mut_env_effects(const int index,
    std::vector<int>& ret_key, std::vector<double>& ret_value) const;
  void return_all_indices(std::vector<int>& ret_index_genotype,
    std::vector<std::vector<int>>& ret_index_interactions,
    std::vector<int>& ret_index_env_effects) const;
  void set_genotype_mutation(const int index,
    const int mut_index, const double effect);
  void set_interaction_mutation(const int index1, const int index2,
    const int mut_index, const double effect);
  void set_env_effects_mutation(const int index,
    const int mut_index, const double effect);
  void set_index_genotype(const std::vector<int>& input);
  void set_index_interactions(const int index, const std::vector<int>& input);
  void set_index_env_effects(const std::vector<int>& input);
};

#endif
