#include <iostream>
#include <fstream>
#include <cmath>
#include "population.hpp"
#include "parameters.hpp"

int main(){
  Parameters para;
  para.pop_size = 1000;

  para.num_g = 2;
  para.alpha = 0.2;

  para.gamma_g = 0.1;
  para.mu_g = 1e-4;
  para.gamma_b = 0.1;
  para.mu_b = 1e-4;
  para.gamma_c = 0.1;
  para.mu_c = 1e-4;

  para.max_g = 1.0 / para.alpha;

  para.sig = 5.0 / std::sqrt(2.0);
  para.criteria = 0.95;
  para.dev_period = 40;
  para.sensor_period = 10;

  para.sep = 1;

  double high = 5.0;
  double low = 0.0;

  para.sig *= (high - low) / para.max_g;

  std::vector<int> tmp_s = {2, -2};

  std::vector<double> s1_vec;
  std::vector<double> s2_vec;

  for(const auto& i: tmp_s){
    if(i == 0){
      s1_vec.push_back(low);
      s2_vec.push_back(low);
    }else if(i == 1){
      s1_vec.push_back(high);
      s2_vec.push_back(high);
    }else if(i == 2){
      s1_vec.push_back(high);
      s2_vec.push_back(low);
    }else if(i == -2){
      s1_vec.push_back(low);
      s2_vec.push_back(high);
    }
  }

  double s1_env = para.max_g;
  double s2_env = 0.0;

  std::vector<double> env_cues = {s1_env, s2_env};
  std::vector<std::vector<double>> s_vecs = {s1_vec, s2_vec};

  Population pop(para, env_cues, s_vecs);

  std::ofstream ofs("regi_fit.txt");

  for(int i = 0; i <= 1e+6; i++){
    if(i % 2000 == 0){
      pop.output_ind();
    }

    pop.one_generation();

    double fit1, fit2;
    pop.return_mean_fitness(fit1, fit2);

    if(i % 100 == 0){
      ofs << i << "\t" << fit1 << "\t" << fit2 << std::endl;
    }
  }

  return(0);
}
