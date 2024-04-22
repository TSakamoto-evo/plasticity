#include "population.hpp"

Population::Population(const Parameters input_para,
  const std::vector<double>& input_env_cues,
  const std::vector<std::vector<double>>& input_s_vecs){

  para = input_para;
  env_cues = input_env_cues;
  s_vecs = input_s_vecs;

  std::vector<double> input_genotype;
  for(int i = 0; i < para.num_g; i++){
    input_genotype.push_back(para.max_g / 2.0);
  }

  std::vector<std::vector<double>> input_interactions;
  for(int i = 0; i < para.num_g; i++){
    std::vector<double> tmp;
    for(int j = 0; j < para.num_g; j++){
      tmp.push_back(0.0);
    }
    input_interactions.push_back(tmp);
  }

  std::vector<double> input_env_effects;
  for(int i = 0; i < para.num_g; i++){
    input_env_effects.push_back(0.0);
  }

  pop.clear();
  pop.emplace_back(para.pop_size, input_genotype, input_interactions, input_env_effects);
  max_index = 0;
  mode = 1;

  next_gen.clear();
  for(int i = 0; i < para.pop_size; i++){
    next_gen.emplace_back(0, input_genotype, input_interactions, input_env_effects);
  }

  std::random_device seed;
  std::mt19937 tmp_mt(seed());
  mt = tmp_mt;

  gen = 0;

  std::poisson_distribution<> tmp_mut_num(para.pop_size * (para.num_g * para.mu_g +
    para.num_g * para.num_g * para.mu_b + para.num_g * para.mu_c));
  std::vector<double> prob = {para.num_g * para.mu_g,
    para.num_g * para.num_g * para.mu_b, para.num_g * para.mu_c};
  std::discrete_distribution<> tmp_choose_type(prob.begin(), prob.end());
  std::uniform_int_distribution<> tmp_choose_locus(0, para.num_g - 1);
  std::uniform_int_distribution<> tmp_choose_uni_ind(0, para.pop_size - 1);
  std::uniform_real_distribution<> tmp_d_gamma_g(-para.gamma_g, para.gamma_g);
  std::uniform_real_distribution<> tmp_d_gamma_b(-para.gamma_b, para.gamma_b);
  std::uniform_real_distribution<> tmp_d_gamma_c(-para.gamma_c, para.gamma_c);

  mut_num = tmp_mut_num;
  choose_type = tmp_choose_type;
  choose_uni_ind = tmp_choose_uni_ind;
  choose_locus = tmp_choose_locus;
  d_gamma_g = tmp_d_gamma_g;
  d_gamma_b = tmp_d_gamma_b;
  d_gamma_c = tmp_d_gamma_c;

  bit_size = 0;

  fitness1.clear();
  fitness2.clear();
  fitness1.reserve(para.pop_size);
  fitness2.reserve(para.pop_size);

  fluct = 0;
}

void Population::development(){
  double deno = 2.0 * para.sig * para.sig;
  mean_fitness1 = 0.0;
  mean_fitness2 = 0.0;

  fitness1.clear();
  fitness2.clear();

  for(int i = 0; i <= max_index; i++){
    fitness1.push_back(pop[i].ret_fitness(para.alpha, para.dev_period, env_cues[0],
      para.sensor_period, s_vecs[0], deno));
    mean_fitness1 += fitness1[i];
  }

  for(int i = 0; i <= max_index; i++){
    fitness2.push_back(pop[i].ret_fitness(para.alpha, para.dev_period, env_cues[1],
      para.sensor_period, s_vecs[1], deno));
    mean_fitness2 += fitness2[i];
  }

  mean_fitness1 /= para.pop_size;
  mean_fitness2 /= para.pop_size;
}

bool Population::return_bit(){
  if(bit_size == 0){
    std::uint32_t result = mt();
    std::bitset<32> tmp(result);
    bs = tmp;
    bit_size += 32;
  }

  bool ret = bs[0];
  bit_size--;
  bs = (bs >> 1);

  return(ret);
}

void Population::selection(){
  std::vector<double> fitness;
  if((gen / para.sep) % 2 == 0){
    fitness = fitness1;
  }else{
    fitness = fitness2;
  }

  std::discrete_distribution<> choose_ind(fitness.begin(), fitness.end());

  std::vector<double> geno1, geno2;
  std::vector<std::vector<double>> inte1, inte2;
  std::vector<double> env1, env2;

  max_index = -1;

  for(int i = 0; i < para.pop_size; i++){
    int parent1 = choose_ind(mt);
    int parent2 = choose_ind(mt);

    pop[parent1].return_genotype(geno1);
    pop[parent1].return_interactions(inte1);
    pop[parent1].return_env_effects(env1);

    pop[parent2].return_genotype(geno2);
    pop[parent2].return_interactions(inte2);
    pop[parent2].return_env_effects(env2);

    for(int j = 0; j < para.num_g; j++){
      if(return_bit()){
        geno1[j] = geno2[j];
      }
      if(return_bit()){
        env1[j] = env2[j];
      }

      for(int k = 0; k < para.num_g; k++){
        if(return_bit()){
          inte1[j][k] = inte2[j][k];
        }
      }
    }

    if(mode == 1 || gen % 50 == 0){
      int success = -1;
      // check if the identical genotype already exists in the nextgen pop
      for(int j = 0; j <= max_index; j++){
        if(next_gen[j].check_identical(geno1, inte1, env1)){
          success = j;
          break;
        }
      }

      if(success >= 0){
        next_gen[success].add_ind();
      }else{
        max_index++;
        next_gen[max_index].set_genotype(1, geno1, inte1, env1);
      }
    }else{
      max_index++;
      next_gen[max_index].set_genotype(1, geno1, inte1, env1);
    }
  }

  if(max_index < 400){
    mode = 1;
  }else{
    mode = 0;
  }

  for(int i = max_index + 1; i < para.pop_size; i++){
    next_gen[i].set_ind_num(0);
  }

  pop = next_gen;
}

void Population::gen_mutation(){
  std::vector<int> ind_num_dist;
  for(int i = 0; i <= max_index; i++){
    ind_num_dist.push_back(pop[i].ret_ind_num());
  }

  int num_mut = mut_num(mt);

  for(int i = 0; i < num_mut; i++){
    int ind = choose_uni_ind(mt);

    int index = 0;
    while(ind >= ind_num_dist.at(index)){
      ind -= ind_num_dist.at(index);
      index++;
    }

    if(ind_num_dist.at(index) > 1){
      ind_num_dist.at(index)--;
      ind_num_dist.push_back(1);

      pop.at(index).reduce_ind();
      max_index++;
      pop.at(max_index) = pop.at(index);
      pop.at(max_index).set_ind_num(1);
      index = max_index;
    }

    int type = choose_type(mt);

    if(type == 0){
      int locus = choose_locus(mt);
      pop.at(index).add_mut(locus, d_gamma_g(mt), para.max_g);
    }else if(type == 1){
      int locus1 = choose_locus(mt);
      int locus2 = choose_locus(mt);
      int locus = para.num_g * locus1 + locus2 + para.num_g;

      pop.at(index).add_mut(locus, d_gamma_b(mt), para.max_g);
    }else{
      int locus1 = choose_locus(mt);
      int locus = locus1 + para.num_g * (para.num_g + 1);
      pop.at(index).add_mut(locus, d_gamma_c(mt), para.max_g);
    }
  }
}

void Population::one_generation(){
  development();
  selection();
  gen_mutation();

  gen++;
}

void Population::make_output_file(const int var_run_index,
  const int para_index, const int rep, const int last_gen, const int suf_time,
  const int max_time) const{

  std::ofstream ofs("regi_state" + std::to_string(var_run_index) + ".txt");
  ofs << var_run_index << "\n" << para_index << "\n" << rep << "\n" << gen <<
    "\n" << last_gen << "\n" << suf_time << "\n" << max_time << "\n" << fluct << std::endl;

  ofs << para.pop_size << "\n" << para.num_g << "\n" << para.alpha << "\n" <<
    para.gamma_g << "\n" << para.mu_g << "\n" << para.gamma_b << "\n" << para.mu_b <<
    "\n" << para.gamma_c << "\n" << para.mu_c << "\n" << para.max_g << "\n" <<
    para.sig << "\n" << para.criteria << "\n" << para.dev_period << "\n" <<
    para.sep << "\n" << para.sensor_period << "\n" << para.num_trial << std::endl;

  ofs << 0 << "\t" << env_cues.at(0) << "\t" << env_cues.at(1) << std::endl;
  ofs << 1;
  for(const auto i: s_vecs.at(0)){
    ofs << "\t" << i;
  }
  ofs << std::endl;

  ofs << 2;
  for(const auto i: s_vecs.at(1)){
    ofs << "\t" << i;
  }
  ofs << std::endl;

  std::vector<double> regi_g;
  std::vector<std::vector<double>> regi_b;
  std::vector<double> regi_c;

  for(int i = 0; i <= max_index; i++){
    pop.at(i).return_genotype(regi_g);
    pop.at(i).return_interactions(regi_b);
    pop.at(i).return_env_effects(regi_c);

    ofs << 3 << "\t" << i << "\t" << 0 << "\t" << pop.at(i).ret_ind_num() << std::endl;

    ofs << 3 << "\t" << i << "\t" << 1;
    for(const auto& j: regi_g){
      ofs << "\t" << j;
    }
    ofs << std::endl;

    for(const auto& j: regi_b){
      ofs << 3 << "\t" << i << "\t" << 2;
      for(const auto& k: j){
        ofs << "\t" << k;
      }
      ofs << std::endl;
    }

    ofs << 3 << "\t" << i << "\t" << 3;
    for(const auto& j: regi_c){
      ofs << "\t" << j;
    }
    ofs << std::endl;
  }

  ofs.close();
}

void Population::check_attractor(std::vector<int>& ret){
  int count1 = 0;
  int count2 = 0;
  int count3 = 0;

  std::uniform_real_distribution<> set_ini_g(0.0, para.max_g);

  std::vector<int> ind_num_dist;
  for(int i = 0; i <= max_index; i++){
    ind_num_dist.push_back(pop[i].ret_ind_num());
  }
  std::discrete_distribution<> choose_ind(ind_num_dist.begin(), ind_num_dist.end());

  for(int i = 0; i < para.num_trial; i++){
    std::vector<double> ini_g;

    for(int j = 0; j < para.num_g; j++){
      ini_g.push_back(set_ini_g(mt));
    }

    Genotype test_ind = pop.at(choose_ind(mt));

    std::vector<std::vector<double>> phenotypes;
    test_ind.return_final(ini_g, phenotypes, para.alpha, para.dev_period);

    double tmp_fitness1 = 1.0;
    double tmp_fitness2 = 1.0;
    int adult_step = static_cast<int>(phenotypes.size());

    for(const auto& phenotype: phenotypes){
      double distance1 = 0.0;
      double distance2 = 0.0;

      for(int j = 0; j < para.num_g; j++){
        distance1 += (phenotype.at(j) - s_vecs.at(0).at(j)) *
          (phenotype.at(j) - s_vecs.at(0).at(j));
        distance2 += (phenotype.at(j) - s_vecs.at(1).at(j)) *
          (phenotype.at(j) - s_vecs.at(1).at(j));
      }

      tmp_fitness1 *= std::exp(-distance1 / 2.0 / para.sig / para.sig / adult_step);
      tmp_fitness2 *= std::exp(-distance2 / 2.0 / para.sig / para.sig / adult_step);
    }

    if(tmp_fitness1 > para.criteria){
      count1++;
    }else if(tmp_fitness2 > para.criteria){
      count2++;
    }else{
      count3++;
    }
  }

  std::vector<int> tmp = {count1, count2, count3};
  ret = tmp;
}

void Population::check_env_grad(const int num_env, std::vector<int>& ret){

  int count1 = 0;
  int count2 = 0;
  int count3 = 0;

  for(int i = 0; i <= num_env; i++){
    std::vector<std::vector<double>> phenotypes;

    for(int j = 0; j <= max_index; j++){
      pop[j].return_final_grad(1.0 * i / num_env * para.max_g, phenotypes, para.alpha,
        para.dev_period, para.sensor_period);

      double tmp_fitness1 = 1.0;
      double tmp_fitness2 = 1.0;
      int adult_step = static_cast<int>(phenotypes.size());

      for(const auto& phenotype: phenotypes){
        double distance1 = 0.0;
        double distance2 = 0.0;

        for(int j = 0; j < para.num_g; j++){
          distance1 += (phenotype.at(j) - s_vecs.at(0).at(j)) *
            (phenotype.at(j) - s_vecs.at(0).at(j));
          distance2 += (phenotype.at(j) - s_vecs.at(1).at(j)) *
            (phenotype.at(j) - s_vecs.at(1).at(j));
        }

        tmp_fitness1 *= std::exp(-distance1 / 2.0 / para.sig / para.sig / adult_step);
        tmp_fitness2 *= std::exp(-distance2 / 2.0 / para.sig / para.sig / adult_step);
      }

      if(tmp_fitness1 > para.criteria){
        count1 += pop[j].ret_ind_num();
      }else if(tmp_fitness2 > para.criteria){
        count2 += pop[j].ret_ind_num();
      }else{
        count3 += pop[j].ret_ind_num();
      }
    }
  }

  std::vector<int> tmp = {count1, count2, count3};
  ret = tmp;
}

void Population::check_fluct(const int wait){
  bool ret = pop.at(0).fluct(wait, env_cues, para.alpha, para.dev_period, 
    para.sensor_period);

  if(ret){
    fluct++;
  }
}
