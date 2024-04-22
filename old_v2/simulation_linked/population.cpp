#include "population.hpp"

Population::Population(const Parameters input_para):
  mutation(input_para.num_g){
  para = input_para;

  std::vector<double> tmp_geno;

  tmp_geno.push_back(para.max_g / 2.0);
  for(int i = 0; i < para.num_g; i++){
    tmp_geno.push_back(0.0);
  }
  tmp_geno.push_back(0.0);

  for(int i = 0; i < para.num_g; i++){
    mutation.add_mutation(i, tmp_geno);
  }

  std::vector<double> input_genotype;
  std::vector<std::vector<double>> input_interactions;
  std::vector<double> input_env_effects;

  std::vector<int> input_mut_allele;

  for(int i = 0; i < para.num_g; i++){
    std::vector<double> tmp_genotype = mutation.return_mut_effect(i, 1);
    input_genotype.push_back(tmp_genotype.at(0));

    std::vector<double> tmp1;
    for(int j = 1; j <= para.num_g; j++){
      tmp1.push_back(tmp_genotype.at(j));
    }
    input_interactions.push_back(tmp1);
    input_env_effects.push_back(tmp_genotype.at(para.num_g + 1));
    input_mut_allele.push_back(1);
  }

  pop.clear();
  pop.emplace_back(input_genotype, input_interactions, input_env_effects,
    input_mut_allele, para.pop_size);

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
  std::uniform_real_distribution<> tmp_d_gamma_g(-para.gamma_g, para.gamma_g);
  std::uniform_real_distribution<> tmp_d_gamma_b(-para.gamma_b, para.gamma_b);
  std::uniform_real_distribution<> tmp_d_gamma_c(-para.gamma_c, para.gamma_c);

  mut_num = tmp_mut_num;
  choose_type = tmp_choose_type;
  choose_locus = tmp_choose_locus;
  d_gamma_g = tmp_d_gamma_g;
  d_gamma_b = tmp_d_gamma_b;
  d_gamma_c = tmp_d_gamma_c;
}

void Population::development(const double env_cue, const std::vector<double> s_vec){
  for(auto& i: pop){
    i.development(para.alpha, para.dev_period, env_cue, para.sensor_period);
    i.set_fitness(s_vec, para.sig);
  }
}

void Population::selection(){
  std::vector<int> locus_index;
  std::vector<int> locus_allele_num;
  std::vector<std::vector<int>> locus_allele_index;
  std::vector<std::vector<std::vector<double>>> locus_allele_effect;
  std::vector<std::vector<double>> locus_freq;
  int total_geno = 1;

  for(int i = 0; i < para.num_g; i++){
    int allele_num = mutation.allele_num(i);
    if(allele_num > 1){
      std::vector<std::vector<double>> tmp0;

      locus_index.push_back(i);
      locus_allele_num.push_back(allele_num);
      locus_allele_index.push_back(mutation.return_allele_index(i, tmp0));
      locus_allele_effect.push_back(tmp0);

      std::vector<double> tmp1(allele_num);
      locus_freq.push_back(tmp1);

      total_geno *= allele_num;
    }
  }

  for(const auto& i: pop){
    for(int j = 0; j < static_cast<int>(locus_index.size()); j++){
      int allele_index;
      double freq = i.return_geno_value(locus_index.at(j), allele_index);

      int k = 0;
      while(locus_allele_index.at(j).at(k) != allele_index){
        k++;
      }
      locus_freq.at(j).at(k) += freq;
    }
  }

  for(auto& i: locus_freq){
    double sum = 0.0;
    for(const auto& j: i){
      sum += j;
    }

    for(auto& j: i){
      j /= sum;
    }
  }

  std::vector<std::vector<int>> index_list(total_geno);
  std::vector<std::vector<int>> index_index_list(total_geno);
  std::vector<std::vector<std::vector<double>>> effect_list(total_geno);
  std::vector<double> freq_list(total_geno);

  for(int i = 0; i < total_geno; i++){
    int div = total_geno;

    std::vector<int> tmp;
    double freq = 1.0;

    for(int j = 0; j < static_cast<int>(locus_allele_num.size()); j++){
      div /= locus_allele_num.at(j);
      int allele_index =
        locus_allele_index.at(j).at((i / div) % locus_allele_num.at(j));
      index_list.at(i).push_back(allele_index);

      index_index_list.at(i).push_back((i / div) % locus_allele_num.at(j));

      std::vector<double> allele_effect =
        locus_allele_effect.at(j).at((i / div) % locus_allele_num.at(j));
      effect_list.at(i).push_back(allele_effect);

      freq *= locus_freq.at(j).at((i / div) % locus_allele_num.at(j));
    }

    freq_list.at(i) = freq;
  }

  {
    double sum = 0.0;
    for(const auto& i: freq_list){
      sum += i;
    }

    for(auto& i: freq_list){
      i /= sum;
    }
  }

  std::vector<int> dist(total_geno);
  std::vector<int> indices(total_geno);
  std::iota(indices.begin(), indices.end(), 0);
  std::shuffle(indices.begin(), indices.end(), mt);

  int remaining_pop = para.pop_size;

  for(int i = 0; i < total_geno - 1; i++){
    double sum_freq = 0.0;
    for(int j = i; j < total_geno; j++){
      sum_freq += freq_list.at(indices.at(j));
    }

    std::binomial_distribution<>
      det_num(remaining_pop, freq_list.at(indices.at(i)) / sum_freq);
    int tmp_num = det_num(mt);
    dist.at(indices.at(i)) = tmp_num;
    remaining_pop -= tmp_num;
  }
  dist.at(indices.at(total_geno - 1)) = remaining_pop;

  std::vector<Genotype> next_gen;
  Genotype base = pop.at(0);

  for(int i = 0; i < total_geno; i++){
    if(dist.at(i) > 0){
      base.set_mutation(locus_index, index_list.at(i), effect_list.at(i),
        dist.at(i));
      next_gen.push_back(base);
    }
  }
  pop = next_gen;

  std::vector<std::vector<int>> exist(locus_allele_index);
  for(auto& i: exist){
    for(auto& j: i){
      j = 0;
    }
  }

  for(int i = 0; i < total_geno; i++){
    if(dist.at(i) > 0){
      for(int j = 0; j < static_cast<int>(locus_allele_num.size()); j++){
        exist.at(j).at(index_index_list.at(i).at(j)) = 1;
      }
    }
  }

  for(int i = 0; i < static_cast<int>(exist.size()); i++){
    for(int j = 0; j < static_cast<int>(exist.at(i).size()); j++){
      if(exist.at(i).at(j) == 0){
        mutation.erase(locus_index.at(i), locus_allele_index.at(i).at(j));
      }
    }
  }
}

void Population::gen_mutation(){
  int num_mut = mut_num(mt);

  for(int i = 0; i < num_mut; i++){
    std::vector<int> pop_dist;
    for(const auto& j: pop){
      pop_dist.push_back(j.return_size());
    }
    std::discrete_distribution<> choose_ind(pop_dist.begin(), pop_dist.end());

    int geno_index = choose_ind(mt);
    int type = choose_type(mt);

    if(type == 0){
      int locus = choose_locus(mt);
      Genotype add = pop.at(geno_index);
      pop.at(geno_index).reduce_ind();

      int allele_index = mutation.return_index(locus) + 1;
      std::vector<double> effect = add.add_mut(locus, d_gamma_g(mt), allele_index, para.max_g);
      mutation.add_mutation(locus, effect);
      pop.push_back(add);
    }else if(type == 1){
      int locus1 = choose_locus(mt);
      int locus2 = choose_locus(mt);
      int locus = para.num_g * locus1 + locus2 + para.num_g;

      Genotype add = pop.at(geno_index);
      pop.at(geno_index).reduce_ind();

      int allele_index = mutation.return_index(locus1) + 1;
      std::vector<double> effect = add.add_mut(locus, d_gamma_b(mt), allele_index, para.max_g);
      mutation.add_mutation(locus1, effect);
      pop.push_back(add);
    }else{
      int locus1 = choose_locus(mt);
      int locus = locus1 + para.num_g * (para.num_g + 1);

      Genotype add = pop.at(geno_index);
      pop.at(geno_index).reduce_ind();

      int allele_index = mutation.return_index(locus1) + 1;
      std::vector<double> effect = add.add_mut(locus, d_gamma_c(mt), allele_index, para.max_g);
      mutation.add_mutation(locus1, effect);
      pop.push_back(add);
    }
  }
}

void Population::one_generation(const double env_cue, const std::vector<double> s_vec){
  development(env_cue, s_vec);
  selection();
  gen_mutation();
  gen++;
}

int Population::return_total_geno() const{
  int total_geno = 1;

  for(int i = 0; i < para.num_g * (para.num_g + 2); i++){
    int allele_num = mutation.allele_num(i);
    if(allele_num > 1){
      total_geno *= allele_num;
    }
  }

  return(total_geno);
}

std::vector<double> Population::calculate_mean_fitness(
  const std::vector<double>& env_cues,
  const std::vector<std::vector<double>> s_vecs){

  std::vector<double> ret;

  for(int i = 0; i < static_cast<int>(env_cues.size()); i++){
    development(env_cues.at(i), s_vecs.at(i));

    double sum = 0.0;
    for(const auto& j: pop){
      sum += j.return_fit_size();
    }
    ret.push_back(sum / para.pop_size);
  }

  return(ret);
}

void Population::check_attractor(const std::vector<double>& s1_vec,
  const std::vector<double>& s2_vec, std::vector<int>& ret){

  std::vector<int> pop_dist;
  for(const auto& i: pop){
    pop_dist.push_back(i.return_size());
  }
  std::discrete_distribution<> choose_ind(pop_dist.begin(), pop_dist.end());

  int count1 = 0;
  int count2 = 0;
  int count3 = 0;

  std::uniform_real_distribution<> set_ini_g(0.0, para.max_g);

  for(int i = 0; i < para.num_trial; i++){
    std::vector<double> ini_g;

    for(int j = 0; j < para.num_g; j++){
      ini_g.push_back(set_ini_g(mt));
    }

    Genotype test_ind = pop.at(choose_ind(mt));

    std::vector<double> ret;
    test_ind.return_final(ini_g, ret, para.alpha, para.dev_period);

    double distance1 = 0.0;
    double distance2 = 0.0;

    for(int j = 0; j < para.num_g; j++){
      distance1 += (ret.at(j) - s1_vec.at(j)) * (ret.at(j) - s1_vec.at(j));
      distance2 += (ret.at(j) - s2_vec.at(j)) * (ret.at(j) - s2_vec.at(j));
    }

    double fitness1 = std::exp(-distance1 / 2.0 / para.sig / para.sig);
    double fitness2 = std::exp(-distance2 / 2.0 / para.sig / para.sig);

    if(fitness1 > para.criteria){
      count1++;
    }else if(fitness2 > para.criteria){
      count2++;
    }else{
      count3++;
    }
  }

  std::vector<int> tmp = {count1, count2, count3};
  ret = tmp;
}

void Population::check_env_grad(const std::vector<double>& s1_vec,
  const std::vector<double>& s2_vec, const int num_env,
  std::vector<int>& ret){

  std::vector<int> pop_dist;
  for(const auto& i: pop){
    pop_dist.push_back(i.return_size());
  }
  std::discrete_distribution<> choose_ind(pop_dist.begin(), pop_dist.end());
  Genotype test_ind = pop.at(choose_ind(mt));

  int count1 = 0;
  int count2 = 0;
  int count3 = 0;

  for(int i = 0; i <= num_env; i++){
    std::vector<double> ret;
    test_ind.return_final_grad(1.0 * i / num_env * para.max_g, ret, para.alpha,
      para.dev_period, para.sensor_period);

    double distance1 = 0.0;
    double distance2 = 0.0;

    for(int j = 0; j < para.num_g; j++){
      distance1 += (ret.at(j) - s1_vec.at(j)) * (ret.at(j) - s1_vec.at(j));
      distance2 += (ret.at(j) - s2_vec.at(j)) * (ret.at(j) - s2_vec.at(j));
    }

    double fitness1 = std::exp(-distance1 / 2.0 / para.sig / para.sig);
    double fitness2 = std::exp(-distance2 / 2.0 / para.sig / para.sig);

    if(fitness1 > para.criteria){
      count1++;
    }else if(fitness2 > para.criteria){
      count2++;
    }else{
      count3++;
    }
  }

  std::vector<int> tmp = {count1, count2, count3};
  ret = tmp;
}
