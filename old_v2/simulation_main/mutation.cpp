#include "mutation.hpp"

Mutation::Mutation(const int input_num_g){
  num_g = input_num_g;

  std::vector<std::unordered_map<int, double>> tmp_mut_genotype(num_g,
    std::unordered_map<int, double>());
  std::vector<std::vector<std::unordered_map<int, double>>> tmp_mut_interactions(
    num_g, std::vector<std::unordered_map<int, double>>(num_g,
    std::unordered_map<int, double>())
  );
  std::vector<std::unordered_map<int, double>> tmp_mut_env_effects(num_g,
    std::unordered_map<int, double>());

  mut_genotype = tmp_mut_genotype;
  mut_interactions = tmp_mut_interactions;
  mut_env_effects = tmp_mut_env_effects;

  std::vector<int> tmp_index_genotype(num_g);
  std::vector<std::vector<int>> tmp_index_interactions(num_g,
    std::vector<int>(num_g));
  std::vector<int> tmp_index_env_effects(num_g);

  index_genotype = tmp_index_genotype;
  index_interactions = tmp_index_interactions;
  index_env_effects = tmp_index_env_effects;

  num_genotype = tmp_index_genotype;
  num_interactions = tmp_index_interactions;
  num_env_effects = tmp_index_env_effects;
}

void Mutation::add_mutation(const int index, const double effect){
  if(index < num_g){
    index_genotype.at(index)++;
    int mut_index = index_genotype.at(index);

    mut_genotype.at(index)[mut_index] = effect;
    num_genotype.at(index)++;
  }else if(index < num_g * (num_g + 1)){
    int locus1 = (index - num_g) / num_g;
    int locus2 = (index - num_g) % num_g;

    index_interactions.at(locus1).at(locus2)++;
    int mut_index = index_interactions.at(locus1).at(locus2);

    mut_interactions.at(locus1).at(locus2)[mut_index] = effect;
    num_interactions.at(locus1).at(locus2)++;
  }else{
    index_env_effects.at(index - num_g * (num_g + 1))++;
    int mut_index = index_env_effects.at(index - num_g * (num_g + 1));

    mut_env_effects.at(index - num_g * (num_g + 1))[mut_index] = effect;
    num_env_effects.at(index - num_g * (num_g + 1))++;
  }
}

double Mutation::return_mut_effect(const int index, const int key) const{
  if(index < num_g){
    return(mut_genotype.at(index).at(key));
  }else if(index < num_g * (num_g + 1)){
    int locus1 = (index - num_g) / num_g;
    int locus2 = (index - num_g) % num_g;

    return(mut_interactions.at(locus1).at(locus2).at(key));
  }else{
    return(mut_env_effects.at(index - num_g * (num_g + 1)).at(key));
  }
}

int Mutation::allele_num(const int index) const{
  if(index < num_g){
    return(num_genotype.at(index));
  }else if(index < num_g * (num_g + 1)){
    int locus1 = (index - num_g) / num_g;
    int locus2 = (index - num_g) % num_g;

    return(num_interactions.at(locus1).at(locus2));
  }else{
    return(num_env_effects.at(index - num_g * (num_g + 1)));
  }
}

std::vector<int> Mutation::return_allele_index(const int index,
  std::vector<double>& effects) const{

  std::vector<int> ret;
  effects.clear();

  if(index < num_g){
    for(const auto& i: mut_genotype.at(index)){
      ret.push_back(i.first);
      effects.push_back(i.second);
    }
    return(ret);
  }else if(index < num_g * (num_g + 1)){
    int locus1 = (index - num_g) / num_g;
    int locus2 = (index - num_g) % num_g;

    for(const auto& i: mut_interactions.at(locus1).at(locus2)){
      ret.push_back(i.first);
      effects.push_back(i.second);
    }
    return(ret);
  }else{
    for(const auto& i: mut_env_effects.at(index - num_g * (num_g + 1))){
      ret.push_back(i.first);
      effects.push_back(i.second);
    }
    return(ret);
  }
}

void Mutation::erase(const int locus_index, const int allele_index){
  if(locus_index < num_g){
    mut_genotype.at(locus_index).erase(allele_index);
    num_genotype.at(locus_index)--;
  }else if(locus_index < num_g * (num_g + 1)){
    int locus1 = (locus_index - num_g) / num_g;
    int locus2 = (locus_index - num_g) % num_g;

    mut_interactions.at(locus1).at(locus2).erase(allele_index);
    num_interactions.at(locus1).at(locus2)--;
  }else{
    mut_env_effects.at(locus_index - num_g * (num_g + 1)).erase(allele_index);
    num_env_effects.at(locus_index - num_g * (num_g + 1))--;
  }
}

int Mutation::return_index(const int index) const{
  if(index < num_g){
    return(index_genotype.at(index));
  }else if(index < num_g * (num_g + 1)){
    int locus1 = (index - num_g) / num_g;
    int locus2 = (index - num_g) % num_g;

    return(index_interactions.at(locus1).at(locus2));
  }else{
    return(index_env_effects.at(index - num_g * (num_g + 1)));
  }
}

void Mutation::check_status() const{
  for(int i = 0; i < num_g; i++){
    int num = num_genotype.at(i);

    int num2 = mut_genotype.at(i).size();

    if(num != num2){
      std::cout << "g" << i << " error!" << std::endl;
    }
  }

  for(int i = 0; i < num_g; i++){
    for(int j = 0; j < num_g; j++){
      int num = num_interactions.at(i).at(j);

      int num2 = mut_interactions.at(i).at(j).size();

      if(num != num2){
        std::cout << "b" << i << j << " error!" << std::endl;
      }
    }
  }

  for(int i = 0; i < num_g; i++){
    int num = num_env_effects.at(i);

    int num2 = mut_env_effects.at(i).size();

    if(num != num2){
      std::cout << "c" << i << "\t" << num << "\t" << num2 << " error!" << std::endl;
    }
  }
}

void Mutation::ret_mut_genotype(const int index,
  std::vector<int>& ret_key, std::vector<double>& ret_value) const{

  ret_key.clear();
  ret_value.clear();

  for(const auto& i: mut_genotype.at(index)){
    ret_key.push_back(i.first);
    ret_value.push_back(i.second);
  }
}

void Mutation::ret_mut_interaction(const int index1, const int index2,
  std::vector<int>& ret_key, std::vector<double>& ret_value) const{

  ret_key.clear();
  ret_value.clear();

  for(const auto& i: mut_interactions.at(index1).at(index2)){
    ret_key.push_back(i.first);
    ret_value.push_back(i.second);
  }
}

void Mutation::ret_mut_env_effects(const int index,
  std::vector<int>& ret_key, std::vector<double>& ret_value) const{

  ret_key.clear();
  ret_value.clear();

  for(const auto& i: mut_env_effects.at(index)){
    ret_key.push_back(i.first);
    ret_value.push_back(i.second);
  }
}

void Mutation::return_all_indices(std::vector<int>& ret_index_genotype,
  std::vector<std::vector<int>>& ret_index_interactions,
  std::vector<int>& ret_index_env_effects) const{

  ret_index_genotype = index_genotype;
  ret_index_interactions = index_interactions;
  ret_index_env_effects = index_env_effects;
}

void Mutation::set_genotype_mutation(const int index,
  const int mut_index, const double effect){

  mut_genotype.at(index)[mut_index] = effect;
  num_genotype.at(index)++;
}

void Mutation::set_interaction_mutation(const int index1, const int index2,
  const int mut_index, const double effect){

  mut_interactions.at(index1).at(index2)[mut_index] = effect;
  num_interactions.at(index1).at(index2)++;
}

void Mutation::set_env_effects_mutation(const int index,
  const int mut_index, const double effect){

  mut_env_effects.at(index)[mut_index] = effect;
  num_env_effects.at(index)++;
}

void Mutation::set_index_genotype(const std::vector<int>& input){
  index_genotype = input;
}

void Mutation::set_index_interactions(const int index,
  const std::vector<int>& input){
  index_interactions.at(index) = input;
}

void Mutation::set_index_env_effects(const std::vector<int>& input){
  index_env_effects = input;
}
