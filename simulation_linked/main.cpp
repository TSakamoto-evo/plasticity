#include <iostream>
#include <fstream>
#include "population.hpp"
#include "parameters.hpp"

int main(int argc, char *argv[]){
  int var_pop_size, var_dev_period, var_sep, var_run_index;
  double var_alpha, var_gamma_g, var_gamma_b, var_gamma_c, var_sig, var_criteria;
  int var_s1, var_s2, var_s3, var_s4, var_s5, var_s6, var_s7, var_s8;
  int para_index;

  if(argc == 20){
    sscanf(argv[1], "%d", &var_run_index);
    sscanf(argv[2], "%d", &var_pop_size);
    sscanf(argv[3], "%d", &var_dev_period);
    sscanf(argv[4], "%d", &var_sep);
    sscanf(argv[5], "%lf", &var_alpha);
    sscanf(argv[6], "%lf", &var_gamma_g);
    sscanf(argv[7], "%lf", &var_gamma_b);
    sscanf(argv[8], "%lf", &var_gamma_c);
    sscanf(argv[9], "%lf", &var_sig);
    sscanf(argv[10], "%lf", &var_criteria);
    sscanf(argv[11], "%d", &var_s1);
    sscanf(argv[12], "%d", &var_s2);
    sscanf(argv[13], "%d", &var_s3);
    sscanf(argv[14], "%d", &var_s4);
    sscanf(argv[15], "%d", &var_s5);
    sscanf(argv[16], "%d", &var_s6);
    sscanf(argv[17], "%d", &var_s7);
    sscanf(argv[18], "%d", &var_s8);
    sscanf(argv[19], "%d", &para_index);
  }else{
    std::cout << "Error! The number of parameters is different!" << std::endl;
    return(1);
  }

  Parameters para;
  para.pop_size = var_pop_size;

  para.num_g = 8;
  para.alpha = var_alpha;

  para.gamma_g = var_gamma_g;
  para.mu_g = 1e-7;
  para.gamma_b = var_gamma_b;
  para.mu_b = 1e-7;
  para.gamma_c = var_gamma_c;
  para.mu_c = 1e-7;

  para.max_g = 1.0 / para.alpha;

  para.sig = var_sig;
  para.criteria = var_criteria;
  para.dev_period = var_dev_period;
  para.sensor_period = 10;

  para.sep = var_sep;

  para.num_trial = 1000;

  int suf_time = 1000 * para.pop_size;
  int max_time = 50000 * para.pop_size;

  std::vector<int> tmp_s = {var_s1, var_s2, var_s3, var_s4, var_s5, var_s6,
    var_s7, var_s8};

  std::vector<double> s1_vec;
  std::vector<double> s2_vec;

  for(const auto& i: tmp_s){
    if(i == 0){
      s1_vec.push_back(0);
      s2_vec.push_back(0);
    }else if(i == 1){
      s1_vec.push_back(para.max_g);
      s2_vec.push_back(para.max_g);
    }else if(i == 2){
      s1_vec.push_back(para.max_g);
      s2_vec.push_back(0);
    }else if(i == -2){
      s1_vec.push_back(0);
      s2_vec.push_back(para.max_g);
    }
  }

  double s1_env = para.max_g;
  double s2_env = 0.0;

  std::vector<double> env_cues = {s1_env, s2_env};
  std::vector<std::vector<double>> s_vecs = {s1_vec, s2_vec};

  Population pop(para);
  int last_gen = 0;
  int succeeed_or_not = -1;

  for(int i = 0; i <= max_time + suf_time; i++){
    if((i / para.sep) % 2 == 0){
      pop.one_generation(s1_env, s1_vec);
    }else{
      pop.one_generation(s2_env, s2_vec);
    }

    std::vector<double> fit;
    fit = pop.calculate_mean_fitness(env_cues, s_vecs);

    if(fit.at(0) < para.criteria || fit.at(1) < para.criteria){
      last_gen = i;
    }

    if(i - last_gen == suf_time){
      std::ofstream ofs2("regi_time.txt", std::ios::app);
      ofs2 << var_run_index << "\t" << para_index << "\t1\t" <<
        (last_gen + 1) * 1.0 / para.pop_size << std::endl;
      ofs2.close();

      succeeed_or_not = 1;

      break;
    }

    if(i == max_time + suf_time){
      std::ofstream ofs2("regi_time.txt", std::ios::app);
      ofs2 << var_run_index << "\t" << para_index << "\t0\t" <<
        max_time * 1.0 / para.pop_size << std::endl;
      ofs2.close();

      succeeed_or_not = 0;

      break;
    }
  }

  std::vector<int> ret;
  pop.check_env_grad(s1_vec, s2_vec, 100, ret);
  std::ofstream ofs5("regi_env_grad.txt", std::ios::app);
  ofs5 << var_run_index << "\t" << para_index << "\t" << succeeed_or_not <<
    "\t" << ret.at(0) << "\t" << ret.at(1) << "\t" << ret.at(2) << std::endl;
  ofs5.close();


  return(0);
}
