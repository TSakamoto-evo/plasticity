#include "population.hpp"

Population::Population(const Parameters input_para):
  mutation(input_para.num_g){

  para = input_para;

  for(int i = 0; i < para.num_g; i++){
    mutation.add_mutation(i, para.max_g / 2.0);
  }
  for(int i = 0; i < (para.num_g * para.num_g); i++){
    mutation.add_mutation(i + para.num_g, 0.0);
  }
  for(int i = 0; i < para.num_g; i++){
    mutation.add_mutation(i + para.num_g * (para.num_g + 1), 0.0);
  }

  std::vector<double> input_genotype;
  std::vector<int> input_mut_genotype;
  for(int i = 0; i < para.num_g; i++){
    input_genotype.push_back(mutation.return_mut_effect(i, 1));
    input_mut_genotype.push_back(1);
  }

  std::vector<std::vector<double>> input_interactions;
  std::vector<std::vector<int>> input_mut_interactions;
  for(int i = 0; i < para.num_g; i++){
    std::vector<double> tmp1;
    std::vector<int> tmp2;
    for(int j = 0; j < para.num_g; j++){
      tmp1.push_back(mutation.return_mut_effect(para.num_g + para.num_g * i + j, 1));
      tmp2.push_back(1);
    }

    input_interactions.push_back(tmp1);
    input_mut_interactions.push_back(tmp2);
  }

  std::vector<double> input_env_effects;
  std::vector<int> input_mut_env_effects;
  for(int i = 0; i < para.num_g; i++){
    input_env_effects.push_back(mutation.return_mut_effect(para.num_g * (para.num_g + 1) + i, 1));
    input_mut_env_effects.push_back(1);
  }

  pop.clear();
  pop.emplace_back(input_genotype, input_interactions, input_env_effects,
    input_mut_genotype, input_mut_interactions, input_mut_env_effects,
    para.pop_size);

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

  std::vector<std::vector<int>> tmp_ind_index_list(para.pop_size,
    std::vector<int>(para.num_g * (para.num_g + 2)));
  std::vector<std::vector<double>> tmp_ind_effect_list(para.pop_size,
    std::vector<double>(para.num_g * (para.num_g + 2)));
  std::vector<long long int> tmp_ind_geno_index(para.pop_size);

  ind_index_list = tmp_ind_index_list;
  ind_effect_list = tmp_ind_effect_list;
  ind_geno_index = tmp_ind_geno_index;

  std::vector<int> tmp_ind_rand_ini(para.pop_size);
  std::iota(tmp_ind_rand_ini.begin(), tmp_ind_rand_ini.end(), 0);
  ind_rand_ini = tmp_ind_rand_ini;

  choose_inds.clear();
  for(int i = 0; i < para.pop_size; i++){
    std::uniform_int_distribution<> tmp(0, para.pop_size - i - 1);
    choose_inds.push_back(tmp);
  }
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
  std::vector<std::vector<double>> locus_allele_effect;
  std::vector<std::vector<double>> locus_freq;
  int total_geno = 1;

  for(int i = 0; i < para.num_g * (para.num_g + 2); i++){
    int allele_num = mutation.allele_num(i);
    if(allele_num > 1){
      std::vector<double> tmp0;

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
  std::vector<std::vector<double>> effect_list(total_geno);
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

      double allele_effect =
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

void Population::selection_many_segregation_v2(){
  std::vector<int> locus_index;
  std::vector<int> locus_allele_num;
  std::vector<std::vector<int>> locus_allele_index;
  std::vector<std::vector<double>> locus_allele_effect;
  std::vector<std::vector<double>> locus_freq;

  for(int i = 0; i < para.num_g * (para.num_g + 2); i++){
    int allele_num = mutation.allele_num(i);
    if(allele_num > 1){
      std::vector<double> tmp0;

      locus_index.push_back(i);
      locus_allele_num.push_back(allele_num);
      locus_allele_index.push_back(mutation.return_allele_index(i, tmp0));
      locus_allele_effect.push_back(tmp0);

      std::vector<double> tmp1(allele_num);
      locus_freq.push_back(tmp1);
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

  std::vector<std::vector<int>> exist(locus_allele_index);

  for(int i = 0; i < static_cast<int>(locus_index.size()); i++){
    std::vector<int> ind_rand(ind_rand_ini);

    int max_index = -1;
    int max_num = -1;

    std::vector<int> nums(locus_freq[i].size());
    int remaining_pop = para.pop_size;
    for(int j = 0; j < static_cast<int>(locus_freq[i].size()) - 1; j++){
      double sum_freq = 0.0;
      for(int k = j; k < static_cast<int>(locus_freq[i].size()); k++){
        sum_freq += locus_freq[i][k];
      }

      std::binomial_distribution<> det_num(remaining_pop,
        locus_freq[i][j] / sum_freq);
      int tmp_num = det_num(mt);
      remaining_pop -= tmp_num;
      nums[j] = tmp_num;
      if(tmp_num > 0){
        exist[i][j] = 1;
      }else{
        exist[i][j] = 0;
      }

      if(tmp_num > max_num){
        max_index = j;
        max_num = tmp_num;
      }
    }

    nums[locus_freq[i].size() - 1] = remaining_pop;
    if(remaining_pop > 0){
      exist[i][static_cast<int>(locus_freq[i].size()) - 1] = 1;
    }else{
      exist[i][static_cast<int>(locus_freq[i].size()) - 1] = 0;
    }

    if(remaining_pop > max_num){
      max_index = static_cast<int>(locus_freq[i].size()) - 1;
      max_num = remaining_pop;
    }

    remaining_pop = para.pop_size;

    for(int j = 0; j < static_cast<int>(locus_freq[i].size()); j++){
      if(j != max_index){
        for(int k = 0; k < nums.at(j); k++){
          int tmp_ind_index = choose_inds[para.pop_size - remaining_pop](mt);
          int tmp_ind = ind_rand[tmp_ind_index];
          ind_rand[tmp_ind_index] = ind_rand.back();
          ind_rand.pop_back();
          remaining_pop--;

          ind_index_list[tmp_ind][i] = locus_allele_index[i][j];
          ind_effect_list[tmp_ind][i] = locus_allele_effect[i][j];

          if(i == 0){
            ind_geno_index[tmp_ind] = j;
          }else{
            ind_geno_index[tmp_ind] *= locus_allele_num[i];
            ind_geno_index[tmp_ind] += j;
          }
        }
      }
    }

    for(const auto& j: ind_rand){
      ind_index_list[j][i] = locus_allele_index[i][max_index];
      ind_effect_list[j][i] = locus_allele_effect[i][max_index];

      if(i == 0){
        ind_geno_index[j] = max_index;
      }else{
        ind_geno_index[j] *= locus_allele_num[i];
        ind_geno_index[j] += max_index;
      }
    }
  }

  std::vector<Genotype> next_gen;
  Genotype base = pop.at(0);

  std::unordered_map<long long int, int> return_pop_pos;

  for(int i = 0; i < para.pop_size; i++){
    if(return_pop_pos.count(ind_geno_index.at(i)) == 1){
      next_gen.at(return_pop_pos.at(ind_geno_index.at(i))).add_ind();
    }else{
      return_pop_pos.emplace(ind_geno_index.at(i), static_cast<int>(next_gen.size()));
      base.set_mutation(locus_index, ind_index_list.at(i), ind_effect_list.at(i), 1);
      next_gen.push_back(base);
    }
  }

  pop = next_gen;

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
      double effect = add.add_mut(locus, d_gamma_g(mt), allele_index, para.max_g);
      mutation.add_mutation(locus, effect);
      pop.push_back(add);
    }else if(type == 1){
      int locus1 = choose_locus(mt);
      int locus2 = choose_locus(mt);
      int locus = para.num_g * locus1 + locus2 + para.num_g;

      Genotype add = pop.at(geno_index);
      pop.at(geno_index).reduce_ind();

      int allele_index = mutation.return_index(locus) + 1;
      double effect = add.add_mut(locus, d_gamma_b(mt), allele_index, para.max_g);
      mutation.add_mutation(locus, effect);
      pop.push_back(add);
    }else{
      int locus1 = choose_locus(mt);
      int locus = locus1 + para.num_g * (para.num_g + 1);

      Genotype add = pop.at(geno_index);
      pop.at(geno_index).reduce_ind();

      int allele_index = mutation.return_index(locus) + 1;
      double effect = add.add_mut(locus, d_gamma_c(mt), allele_index, para.max_g);
      mutation.add_mutation(locus, effect);
      pop.push_back(add);
    }
  }
}

void Population::one_generation(const double env_cue, const std::vector<double> s_vec,
  const int threshold){

  development(env_cue, s_vec);

  if(return_total_geno() < threshold && return_total_geno() > 0){
    selection();
  }else{
    selection_many_segregation_v2();
  }

  gen_mutation();

  gen++;
}

long long int Population::return_total_geno() const{
  long long int total_geno = 1;
  double log_total_geno = 0;

  for(int i = 0; i < para.num_g * (para.num_g + 2); i++){
    int allele_num = mutation.allele_num(i);
    if(allele_num > 1){
      total_geno *= allele_num;
      log_total_geno += std::log(1.0 * allele_num);
    }
  }

  if(log_total_geno < std::log(1.0 * LLONG_MAX)){
    return(total_geno);
  }else{
    return(-1);
  }
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

void Population::make_output_file(const int var_run_index,
  const int para_index, const int time, const std::vector<double> env_cues,
  const std::vector<std::vector<double>> s_vecs,
  const int last_gen, const int suf_time, const int max_time) const{

  std::ofstream ofs("regi_state" + std::to_string(var_run_index) + ".txt");
  ofs << var_run_index << "\n" << para_index << "\n" << time <<
    "\n" << last_gen << "\n" << suf_time << "\n" << max_time << std::endl;

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

  int pop_size;
  std::vector<int> mut_genotype;
  std::vector<std::vector<int>> mut_interactions;
  std::vector<int> mut_env_effects;

  std::vector<int> keys;
  std::vector<double> values;
  for(int i = 0; i < para.num_g; i++){
    mutation.ret_mut_genotype(i, keys, values);

    for(int k = 0; k < static_cast<int>(keys.size()); k++){
      ofs << 3 << "\t" << i << "\t"<< keys.at(k) << "\t" <<
        values.at(k) << std::endl;
    }
  }

  for(int i = 0; i < para.num_g; i++){
    for(int j = 0; j < para.num_g; j++){
      mutation.ret_mut_interaction(i, j, keys, values);

      for(int k = 0; k < static_cast<int>(keys.size()); k++){
        ofs << 4 << "\t" << i << "\t" << j << "\t"<< keys.at(k) << "\t" <<
          values.at(k) << std::endl;
      }
    }
  }

  for(int i = 0; i < para.num_g; i++){
    mutation.ret_mut_env_effects(i, keys, values);

    for(int k = 0; k < static_cast<int>(keys.size()); k++){
      ofs << 5 << "\t" << i << "\t"<< keys.at(k) << "\t" <<
        values.at(k) << std::endl;
    }
  }

  mutation.return_all_indices(mut_genotype, mut_interactions, mut_env_effects);

  ofs << 6;
  for(const auto& j: mut_genotype){
    ofs << "\t" << j;
  }
  ofs << std::endl;

  for(const auto& j: mut_interactions){
    ofs << 7;
    for(const auto k: j){
      ofs << "\t" << k;
    }
    ofs << std::endl;
  }

  ofs << 8;
  for(const auto& j: mut_env_effects){
    ofs << "\t" << j;
  }
  ofs << std::endl;

  for(int i = 0; i < static_cast<int>(pop.size()); i++){
    pop.at(i).return_all_indices(pop_size, mut_genotype,
      mut_interactions, mut_env_effects);

    ofs << 9 << "\t" << i << "\t" << 1 << "\t" << pop_size << std::endl;
    ofs << 9 << "\t" << i << "\t" << 2;
    for(const auto& j: mut_genotype){
      ofs << "\t" << j;
    }
    ofs << std::endl;

    for(const auto& j: mut_interactions){
      ofs << 9 << "\t" << i << "\t" << 3;
      for(const auto& k: j){
        ofs << "\t" << k;
      }
      ofs << std::endl;
    }

    ofs << 9 << "\t" << i << "\t" << 4;
    for(const auto& j: mut_env_effects){
      ofs << "\t" << j;
    }
    ofs << std::endl;
  }

  ofs.close();
}
