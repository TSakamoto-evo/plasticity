#ifndef POPULATION
#define POPULATION

#include "genotype.hpp"
#include "parameters.hpp"
#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <bitset>

class Population{
private:
  std::vector<Genotype> pop;
  std::vector<Genotype> next_gen;

  Parameters para;
  std::mt19937 mt;
  int gen;

  std::bitset<32> bs;
  int bit_size;

  std::vector<double> fitness1;
  std::vector<double> fitness2;

  std::vector<double> env_cues;
  std::vector<std::vector<double>> s_vecs;

  std::uniform_real_distribution<> uni;
  std::poisson_distribution<> mut_num;
  std::discrete_distribution<> choose_type;
  std::uniform_int_distribution<> choose_locus;
  std::uniform_int_distribution<> choose_uni_ind;
  std::uniform_real_distribution<> d_gamma_g;
  std::uniform_real_distribution<> d_gamma_b;
  std::uniform_real_distribution<> d_gamma_c;

  double mean_fitness1;
  double mean_fitness2;

  int max_index;
  int mode;

  int fluct;

public:
  Population(const std::string file_name, Parameters& ret_para,
    int& var_run_index, int& para_index, int& rep, int& last_gen, int& ini_time,
    int& suf_time, int& max_time);
  void development();
  void selection();
  void gen_mutation();
  void one_generation();
  bool return_bit();
  void return_mean_fitness(double& fit1, double& fit2){
    fit1 = mean_fitness1;
    fit2 = mean_fitness2;
  };
  void make_output_file(const int var_run_index,
    const int para_index, const int rep, const int last_gen, const int suf_time,
    const int max_time) const;
  void check_attractor(std::vector<int>& ret);
  void check_env_grad(const int num_env, std::vector<int>& ret);
  void check_fluct(const int wait);
  int ret_fluct() const { return(fluct); };
};

#endif
