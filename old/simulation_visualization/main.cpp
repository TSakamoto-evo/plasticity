#include <iostream>
#include <fstream>
#include <random>
#include "population.hpp"
#include "parameters.hpp"

int main(){
  Parameters para;
  para.pop_size = 10000;

  para.num_g = 8;
  para.alpha = 0.2;

  para.mu_g = 0.1;
  para.p_g = 1e-7;
  para.mu_b = 0.1;
  para.p_b = 1e-7;
  para.mu_c = 0.1;
  para.p_c = 1e-7;

  para.max_g = 1.0 / para.alpha;

  para.sig = 5.0;
  para.dev_period = 10;

  para.sep = 1000000;

  int max_time = 2e+8;

  Population pop(para);

  double s1_env = 1.0;
  double s2_env = 0.0;

  std::vector<double> s1_vec = {0, 0, 5, 5, 5, 5, 0, 0};
  std::vector<double> s2_vec = {0, 0, 5, 5, 0, 0, 5, 5};

  std::vector<double> env_cues = {s1_env, s2_env};
  std::vector<std::vector<double>> s_vecs = {s1_vec, s2_vec};

  std::ofstream ofs("regi_genotypes.txt");
  std::ofstream ofs2("answers.txt");
  std::ofstream ofs4("fitness.txt");

  ofs2 << "vec1";
  for(const auto& i: s1_vec){
    ofs2 << "\t" << i;
  }
  ofs2 << std::endl;

  ofs2 << "vec2";
  for(const auto& i: s2_vec){
    ofs2 << "\t" << i;
  }
  ofs2 << std::endl;

  for(int i = 0; i <= max_time; i++){
    if(i % (10 * para.pop_size) == 0 || i % (10 * para.pop_size) == para.sep){
      std::vector<double> ret_genotypes;
      std::vector<std::vector<double>> ret_interactions;
      std::vector<double> ret_env_effects;
      std::vector<double> ret_adult1, ret_adult2;

      pop.return_ind(ret_genotypes, ret_interactions, ret_env_effects,
        ret_adult1, ret_adult2, s1_env, s2_env);

      ofs << i;
      for(const auto& j: ret_genotypes){
        ofs << "\t" << j;
      }

      for(const auto& j: ret_interactions){
        for(const auto& k: j){
          ofs << "\t" << k;
        }
      }

      for(const auto& j: ret_env_effects){
        ofs << "\t" << j;
      }

      for(const auto& j: ret_adult1){
        ofs << "\t" << j;
      }
      for(const auto& j: ret_adult2){
        ofs << "\t" << j;
      }

      ofs << std::endl;
    }

    if(i % 10000 == 0){
      std::vector<double> fit;
      fit = pop.calculate_mean_fitness(env_cues, s_vecs);
      ofs4 << i << "\t" << fit.at(0) << "\t" << fit.at(1) << std::endl;
    }

    if((i / para.sep) % 2 == 0){
      pop.one_generation(s1_env, s1_vec);
    }else{
      pop.one_generation(s2_env, s2_vec);
    }
  }

  return(1);
}
