#include <iostream>
#include <fstream>
#include "population.hpp"
#include "parameters.hpp"

int main(int argc, char *argv[]){
  int var_pop_size, var_dev_period, var_sep, var_run_index;
  double var_alpha, var_factor_g, var_factor_b, var_factor_c, var_sig, var_criteria;
  int var_s1, var_s2, var_s3, var_s4, var_s5, var_s6, var_s7, var_s8;
  int para_index, rep;

  if(argc == 21){
    sscanf(argv[1], "%d", &var_run_index);
    sscanf(argv[2], "%d", &var_pop_size);
    sscanf(argv[3], "%d", &var_dev_period);
    sscanf(argv[4], "%d", &var_sep);
    sscanf(argv[5], "%lf", &var_alpha);
    sscanf(argv[6], "%lf", &var_factor_g);
    sscanf(argv[7], "%lf", &var_factor_b);
    sscanf(argv[8], "%lf", &var_factor_c);
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
    sscanf(argv[20], "%d", &rep);
  }else{
    std::cout << "Error! The number of parameters is different!" << std::endl;
    return(1);
  }

  Parameters para;
  para.pop_size = var_pop_size;

  para.num_g = 8;
  para.alpha = var_alpha;

  para.gamma_g = 0.1;
  para.mu_g = 1e-4 * var_factor_g;
  para.gamma_b = 0.1;
  para.mu_b = 1e-4 * var_factor_b;
  para.gamma_c = 0.1;
  para.mu_c = 1e-4 * var_factor_c;

  para.max_g = 1.0 / para.alpha;

  para.sig = var_sig;
  para.criteria = var_criteria;
  para.dev_period = var_dev_period;
  para.sensor_period = 10;

  para.sep = var_sep;
  para.num_trial = 1000;

  int suf_time = 1e+5;
  int max_time = 1e+7;

  double high = 5.0;
  double low = 0.0;

  para.sig *= (high - low) / para.max_g;

  std::vector<int> tmp_s = {var_s1, var_s2, var_s3, var_s4, var_s5, var_s6,
    var_s7, var_s8};

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
  int last_gen = 0;
  int succeeed_or_not = -1;
  int fluct = 0;

  for(int i = 0; i <= max_time + suf_time; i++){
    if(i > 0 && i % (100 * para.pop_size) == 0){
      pop.make_output_file(var_run_index, para_index, rep, last_gen, suf_time,
        max_time);
    }

    pop.one_generation();

    if(i % 1000 == 0){
      pop.check_fluct(100);
    }

    double fit1, fit2;
    pop.return_mean_fitness(fit1, fit2);

    if(fit1 < para.criteria || fit2 < para.criteria){
      last_gen = i;
    }

    if(i - last_gen == suf_time){
      succeeed_or_not = 1;
      fluct = pop.ret_fluct();
      break;
    }

    if(i == max_time + suf_time){
      succeeed_or_not = 0;
      fluct = pop.ret_fluct();
      break;
    }
  }

  std::vector<int> ret, ret2;
  pop.check_attractor(ret);
  pop.check_env_grad(100, ret2);

  std::ofstream ofs("regi_env_grad.txt", std::ios::app);
  ofs << var_run_index << "\t" << para_index << "\t" << rep << "\t" <<
    succeeed_or_not << "\t" << ret2.at(0) << "\t" << ret2.at(1) << "\t" <<
    ret2.at(2) << std::endl;
  ofs.close();

  if(succeeed_or_not == 1){
    std::ofstream ofs2("regi_time.txt", std::ios::app);
    ofs2 << var_run_index << "\t" << para_index << "\t" << rep << "\t1\t" <<
      (last_gen + 1) << std::endl;
    ofs2.close();
  }else{
    std::ofstream ofs2("regi_time.txt", std::ios::app);
    ofs2 << var_run_index << "\t" << para_index << "\t" << rep << "\t0\t" <<
      max_time << std::endl;
    ofs2.close();
  }

  return(0);
}
