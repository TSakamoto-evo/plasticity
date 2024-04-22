#include "genotype.hpp"

Genotype::Genotype(const std::vector<double> input_genotype,
  const std::vector<std::vector<double>> input_interactions,
  const std::vector<double> input_env_effects,
  const std::vector<int> input_mut_genotype,
  const std::vector<std::vector<int>> input_mut_interactions,
  const std::vector<int> input_mut_env_effects,
  const int input_ind_num){

  genotype = input_genotype;
  interactions = input_interactions;
  env_effects = input_env_effects;

  mut_genotype = input_mut_genotype;
  mut_interactions = input_mut_interactions;
  mut_env_effects = input_mut_env_effects;
  ind_num = input_ind_num;

  num_g = static_cast<int>(genotype.size());
}

void Genotype::development(const double alpha, const int round,
  const double env_cue){

  std::vector<double> p_vec(genotype);

  for(int i = 0; i < round; i++){
    std::vector<double> tmp_p_vec(num_g);
    for(int j = 0; j < num_g; j++){
      double effects = 0.0;
      for(int k = 0; k < num_g; k++){
        effects += interactions.at(j).at(k) * p_vec.at(k);
      }
      effects += env_effects.at(j) * env_cue;

      tmp_p_vec.at(j) = p_vec.at(j) + (0.5 + 0.5 * std::tanh(effects)) -
        alpha * p_vec.at(j);
    }
    p_vec = tmp_p_vec;
  }
  phenotype = p_vec;
}

void Genotype::set_fitness(const std::vector<double> s_vec, const double sig){
  double distance2 = 0.0;
  for(int i = 0; i < num_g; i++){
    distance2 += (s_vec.at(i) - phenotype.at(i)) * (s_vec.at(i) - phenotype.at(i));
  }

  fitness = std::exp(-distance2 / 2.0 / sig / sig);
}

double Genotype::return_geno_value(const int index, int& allele_index) const{
  if(index < num_g){
    allele_index = mut_genotype.at(index);
    return(ind_num * fitness);
  }else if(index < num_g * (num_g + 1)){
    int locus1 = (index - num_g) / num_g;
    int locus2 = (index - num_g) % num_g;

    allele_index = mut_interactions.at(locus1).at(locus2);
    return(ind_num * fitness);
  }else{
    allele_index = mut_env_effects.at(index - num_g * (num_g + 1));
    return(ind_num * fitness);
  }
}

void Genotype::set_mutation(const std::vector<int>& locus_index,
  const std::vector<int>& index_list, const std::vector<double>& effect_list,
  const int input_ind_num){

  ind_num = input_ind_num;

  for(int i = 0; i < static_cast<int>(locus_index.size()); i++){
    int index = locus_index.at(i);
    if(index < num_g){
      genotype.at(index) = effect_list.at(i);
      mut_genotype.at(index) = index_list.at(i);
    }else if(index < num_g * (num_g + 1)){
      int locus1 = (index - num_g) / num_g;
      int locus2 = (index - num_g) % num_g;

      interactions.at(locus1).at(locus2) = effect_list.at(i);
      mut_interactions.at(locus1).at(locus2) = index_list.at(i);
    }else{
      env_effects.at(index - num_g * (num_g + 1)) = effect_list.at(i);
      mut_env_effects.at(index - num_g * (num_g + 1)) = index_list.at(i);
    }
  }
}

double Genotype::add_mut(const int index, const double delta,
  const int allele_index, const double max_g){

  ind_num = 1;

  if(index < num_g){
    mut_genotype.at(index) = allele_index;
    genotype.at(index) += delta;
    if(genotype.at(index) < 0.0){
      genotype.at(index) = 0.0;
    }
    if(genotype.at(index) > max_g){
      genotype.at(index) = max_g;
    }

    return(genotype.at(index));
  }else if(index < num_g * (num_g + 1)){
    int locus1 = (index - num_g) / num_g;
    int locus2 = (index - num_g) % num_g;

    mut_interactions.at(locus1).at(locus2) = allele_index;
    interactions.at(locus1).at(locus2) += delta;

    return(interactions.at(locus1).at(locus2));
  }else{
    mut_env_effects.at(index - num_g * (num_g + 1)) = allele_index;
    env_effects.at(index - num_g * (num_g + 1)) += delta;

    return(env_effects.at(index - num_g * (num_g + 1)));
  }
}
