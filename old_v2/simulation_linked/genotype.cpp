#include "genotype.hpp"

Genotype::Genotype(const std::vector<double> input_genotype,
  const std::vector<std::vector<double>> input_interactions,
  const std::vector<double> input_env_effects,
  const std::vector<int> input_mut_allele,
  const int input_ind_num){

  genotype = input_genotype;
  interactions = input_interactions;
  env_effects = input_env_effects;

  mut_allele = input_mut_allele;
  ind_num = input_ind_num;

  num_g = static_cast<int>(genotype.size());
}

void Genotype::development(const double alpha, const int round,
  const double env_cue, const int sensor_period){

  std::vector<double> p_vec(genotype);

  for(int i = 0; i < round; i++){
    std::vector<double> tmp_p_vec(num_g);
    for(int j = 0; j < num_g; j++){
      double effects = 0.0;
      for(int k = 0; k < num_g; k++){
        effects += interactions.at(j).at(k) * p_vec.at(k);
      }

      if(i < sensor_period){
        effects += env_effects.at(j) * env_cue;
      }

      tmp_p_vec.at(j) = p_vec.at(j) + (0.5 + 0.5 * std::tanh(effects)) -
        alpha * p_vec.at(j);
    }
    p_vec = tmp_p_vec;
  }
  phenotype = p_vec;
}

void Genotype::set_fitness(const std::vector<double> s_vec, const double sig){
  double distance = 0.0;
  for(int i = 0; i < num_g; i++){
    distance += (s_vec.at(i) - phenotype.at(i)) * (s_vec.at(i) - phenotype.at(i));
  }

  fitness = std::exp(-distance / 2.0 / sig / sig);
}

double Genotype::return_geno_value(const int index, int& allele_index) const{
  allele_index = mut_allele.at(index);
  return(ind_num * fitness);
}

void Genotype::set_mutation(const std::vector<int>& locus_index,
  const std::vector<int>& index_list,
  const std::vector<std::vector<double>>& effect_list,
  const int input_ind_num){

  ind_num = input_ind_num;

  for(int i = 0; i < static_cast<int>(locus_index.size()); i++){
    int index = locus_index.at(i);

    genotype.at(index) = effect_list.at(i).at(0);

    for(int j = 0; j < num_g; j++){
      interactions.at(index).at(j) = effect_list.at(i).at(j + 1);
    }

    env_effects.at(index) = effect_list.at(i).at(num_g + 1);

    mut_allele.at(index) = index_list.at(i);
  }
}

std::vector<double> Genotype::add_mut(const int index, const double delta,
  const int allele_index, const double max_g){

  ind_num = 1;
  int index_locus;

  if(index < num_g){
    index_locus = index;

    genotype.at(index) += delta;
    if(genotype.at(index) < 0.0){
      genotype.at(index) = 0.0;
    }
    if(genotype.at(index) > max_g){
      genotype.at(index) = max_g;
    }
  }else if(index < num_g * (num_g + 1)){
    int locus1 = (index - num_g) / num_g;
    int locus2 = (index - num_g) % num_g;

    index_locus = locus1;

    interactions.at(locus1).at(locus2) += delta;
  }else{
    index_locus = index - num_g * (num_g + 1);
    env_effects.at(index_locus) += delta;
  }

  mut_allele.at(index_locus) = allele_index;

  std::vector<double> ret;
  ret.push_back(genotype.at(index_locus));
  for(int i = 0; i < num_g; i++){
    ret.push_back(interactions.at(index_locus).at(i));
  }
  ret.push_back(env_effects.at(index_locus));

  return(ret);
}

void Genotype::return_interactions(std::vector<std::vector<double>>& ret) const{
  ret = interactions;
}


void Genotype::return_final(const std::vector<double>& ini_g,
  std::vector<double>& ret, const double alpha, const int round){

  std::vector<double> p_vec(ini_g);

  for(int i = 0; i < round; i++){
    std::vector<double> tmp_p_vec(num_g);
    for(int j = 0; j < num_g; j++){
      double effects = 0.0;
      for(int k = 0; k < num_g; k++){
        effects += interactions.at(j).at(k) * p_vec.at(k);
      }

      tmp_p_vec.at(j) = p_vec.at(j) + (0.5 + 0.5 * std::tanh(effects)) -
        alpha * p_vec.at(j);
    }
    p_vec = tmp_p_vec;
  }
  ret = p_vec;
}

void Genotype::return_final_grad(const double env_cue,
  std::vector<double>& ret, const double alpha, const int round,
  const int sensor_period){

  std::vector<double> p_vec(genotype);

  for(int i = 0; i < round; i++){
    std::vector<double> tmp_p_vec(num_g);
    for(int j = 0; j < num_g; j++){
      double effects = 0.0;
      for(int k = 0; k < num_g; k++){
        effects += interactions.at(j).at(k) * p_vec.at(k);
      }

      if(i < sensor_period){
        effects += env_effects.at(j) * env_cue;
      }

      tmp_p_vec.at(j) = p_vec.at(j) + (0.5 + 0.5 * std::tanh(effects)) -
        alpha * p_vec.at(j);
    }
    p_vec = tmp_p_vec;
  }
  ret = p_vec;
}
