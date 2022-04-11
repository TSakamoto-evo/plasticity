#ifndef MUTATION
#define MUTATION

#include <vector>
#include <unordered_map>
#include <iostream>

class Mutation{
private:
  std::vector<std::unordered_map<int, std::vector<double>>> mut_allele;
  std::vector<int> index_allele;
  std::vector<int> num_allele;
  int num_g;

public:
  Mutation(const int input_num_g);
  void add_mutation(const int index, const std::vector<double> effect);
  std::vector<double> return_mut_effect(const int index, const int key) const;
  int allele_num(const int index) const;
  std::vector<int> return_allele_index(const int index,
    std::vector<std::vector<double>>& effects) const;
  void erase(const int locus_index, const int allele_index);
  int return_index(const int index) const;
  void check_status() const;
};

#endif
