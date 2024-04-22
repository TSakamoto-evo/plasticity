#include <iostream>
#include <fstream>
#include "population_restart.hpp"
#include "parameters.hpp"

int main(int argc, char *argv[]){
  std::string file_name = "";

  if(argc == 2){
    file_name = argv[1];
  }else{
    std::cout << "Error! The number of parameters is different!" << std::endl;
    return(1);
  }

  Parameters para;

  int ini_time = 0;
  int suf_time = 0;
  int max_time = 0;

  int var_run_index = 0;
  int para_index = 0;
  int rep = 0;

  int last_gen = 0;
  int succeeed_or_not = -1;
  int fluct = 0;

  Population pop(file_name, para, var_run_index, para_index, rep, last_gen,
    ini_time, suf_time, max_time);

  for(int i = ini_time; i <= max_time + suf_time; i++){
    if(i > ini_time && i % (100 * para.pop_size) == 0){
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

  std::ofstream ofs5("regi_env_grad.txt", std::ios::app);
  ofs5 << var_run_index << "\t" << para_index << "\t" << rep << "\t" <<
    succeeed_or_not << "\t" << ret2.at(0) << "\t" << ret2.at(1) << "\t" <<
    ret2.at(2) << std::endl;
  ofs5.close();

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
