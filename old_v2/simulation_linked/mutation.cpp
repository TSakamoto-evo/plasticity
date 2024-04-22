#include "mutation.hpp"

Mutation::Mutation(const int input_num_g){
  num_g = input_num_g;

  std::vector<std::unordered_map<int, std::vector<double>>> tmp_mut_allele(num_g,
    std::unordered_map<int, std::vector<double>>());

  mut_allele = tmp_mut_allele;

  std::vector<int> tmp_index_allele(num_g);
  index_allele = tmp_index_allele;
  num_allele = tmp_index_allele;
}

void Mutation::add_mutation(const int index, const std::vector<double> effect){
    index_allele.at(index)++;
    int mut_index = index_allele.at(index);

    mut_allele.at(index)[mut_index] = effect;
    num_allele.at(index)++;
}

std::vector<double> Mutation::return_mut_effect(const int index, const int key) const{
  return(mut_allele.at(index).at(key));
}

int Mutation::allele_num(const int index) const{
    return(num_allele.at(index));
}

std::vector<int> Mutation::return_allele_index(const int index,
  std::vector<std::vector<double>>& effects) const{

  std::vector<int> ret;
  effects.clear();

  for(const auto& i: mut_allele.at(index)){
    ret.push_back(i.first);
    effects.push_back(i.second);
  }
  return(ret);
}

void Mutation::erase(const int locus_index, const int allele_index){
    mut_allele.at(locus_index).erase(allele_index);
    num_allele.at(locus_index)--;
}

int Mutation::return_index(const int index) const{
  return(index_allele.at(index));
}

void Mutation::check_status() const{
  for(int i = 0; i < num_g; i++){
    int num = num_allele.at(i);

    int num2 = mut_allele.at(i).size();

    if(num != num2){
      std::cout << "g" << i << " error!" << std::endl;
    }
  }
}
