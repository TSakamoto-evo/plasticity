#include "genotype.hpp"

Genotype::Genotype(const int input_ind_num,
  const std::vector<double>& input_genotype,
  const std::vector<std::vector<double>>& input_interactions,
  const std::vector<double>& input_env_effects){

  ind_num = input_ind_num;

  genotype = input_genotype;
  interactions = input_interactions;
  env_effects = input_env_effects;

  num_g = static_cast<int>(genotype.size());

  adult_step = 20;
}

void Genotype::set_genotype(const int input_ind_num,
  const std::vector<double>& input_genotype,
  const std::vector<std::vector<double>>& input_interactions,
  const std::vector<double>& input_env_effects){

  ind_num = input_ind_num;

  genotype = input_genotype;
  interactions = input_interactions;
  env_effects = input_env_effects;

  num_g = static_cast<int>(genotype.size());
}

void Genotype::development(const double alpha, const int round,
  const double env_cue, const int sensor_period){

  phenotypes.clear();

  std::vector<double> p_vec(genotype);
  std::vector<double> tmp_p_vec(num_g);

  for(int i = 0; i < round; ++i){
    for(int j = 0; j < num_g; ++j){
      double effects = 0.0;
      for(int k = 0; k < num_g; ++k){
        effects += interactions[j][k] * p_vec[k];
      }

      if(i < sensor_period){
        effects += env_effects[j] * env_cue;
      }

      tmp_p_vec[j] = (1.0 - alpha) * p_vec[j] + (0.5 + 0.5 * std::tanh(effects));
    }
    p_vec = tmp_p_vec;

    if(i >= round - adult_step){
      phenotypes.push_back(p_vec);
    }
  }
}

double Genotype::ret_fitness(const double alpha, const int round, const double env_cue,
  const int sensor_period, const std::vector<double>& s_vec, const double deno){

  development(alpha, round, env_cue, sensor_period);
  double fitness = 1.0;

  for(const auto& phenotype: phenotypes){
    double distance_sq = 0.0;
    for(int i = 0; i < num_g; i++){
      distance_sq += (s_vec[i] - phenotype[i]) * (s_vec[i] - phenotype[i]);
    }
    fitness *= std::exp(-distance_sq / deno / adult_step);
  }
  
  return(fitness * ind_num);
}


void Genotype::add_mut(const int index, const double delta, const double max_g){
  if(index < num_g){
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

    interactions.at(locus1).at(locus2) += delta;
  }else{
    env_effects.at(index - num_g * (num_g + 1)) += delta;
  }
}

void Genotype::return_final(const std::vector<double>& ini_g,
  std::vector<std::vector<double>>& ret, const double alpha, 
  const int round) const{

  ret.clear();

  std::vector<double> p_vec(ini_g);
  std::vector<double> tmp_p_vec(num_g);

  for(int i = 0; i < round; i++){
    for(int j = 0; j < num_g; j++){
      double effects = 0.0;
      for(int k = 0; k < num_g; k++){
        effects += interactions.at(j).at(k) * p_vec.at(k);
      }
      tmp_p_vec.at(j) = p_vec.at(j) + (0.5 + 0.5 * std::tanh(effects)) -
        alpha * p_vec.at(j);
    }
    p_vec = tmp_p_vec;

    if(i >= round - adult_step){
      ret.push_back(p_vec);
    }
  }
}

void Genotype::return_final_grad(const double env_cue,
  std::vector<std::vector<double>>& ret, const double alpha, 
  const int round, const int sensor_period) const{

  ret.clear();

  std::vector<double> p_vec(genotype);
  std::vector<double> tmp_p_vec(num_g);

  for(int i = 0; i < round; i++){
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

    if(i >= round - adult_step){
      ret.push_back(p_vec);
    }
  }
}

bool Genotype::check_identical(const std::vector<double>& input_genotype,
  const std::vector<std::vector<double>>& input_interactions,
  const std::vector<double>& input_env_effects) const{

  for(int i = 0; i < num_g; i++){
    if(genotype[i] != input_genotype[i]){
      return(0);
    }
  }

  for(int i = 0; i < num_g; i++){
    for(int j = 0; j < num_g; j++){
      if(interactions[i][j] != input_interactions[i][j]){
        return(0);
      }
    }
  }

  for(int i = 0; i < num_g; i++){
    if(env_effects[i] != input_env_effects[i]){
      return(0);
    }
  }

  return(1);
}

bool Genotype::fluct(const int wait, const std::vector<double>& env_cues,
    const double alpha, const int round, const int sensor_period) const{

  bool fluct = 0;

  for(const auto& env_cue: env_cues){
    std::vector<double> p_vec(genotype);
    std::vector<double> tmp_p_vec(num_g);

    for(int i = 0; i < wait; i++){
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

    std::vector<double> min_vals(p_vec);
    std::vector<double> max_vals(p_vec);

    for(int i = 0; i < round; i++){
      for(int j = 0; j < num_g; j++){
        double effects = 0.0;
        for(int k = 0; k < num_g; k++){
          effects += interactions.at(j).at(k) * p_vec.at(k);
        }

        tmp_p_vec.at(j) = p_vec.at(j) + (0.5 + 0.5 * std::tanh(effects)) -
          alpha * p_vec.at(j);
      }
      p_vec = tmp_p_vec;

      for(int j = 0; j < num_g; j++){
        if(p_vec.at(j) > max_vals.at(j)){
          max_vals.at(j) = p_vec.at(j);
        }else if(p_vec.at(j) < min_vals.at(j)){
          min_vals.at(j) = p_vec.at(j);
        }
      }
    }

    std::vector<double> mid_vals(num_g);
    std::vector<bool> ample(num_g);
    std::vector<bool> sign(num_g);
    for(int i = 0; i < num_g; i++){
      mid_vals.at(i) = (max_vals.at(i) + min_vals.at(i)) / 2.0;

      if(max_vals.at(i) - min_vals.at(i) > 0.1){
        ample.at(i) = 1;
      }

      if(p_vec.at(i) >= mid_vals.at(i)){
        sign.at(i) = 1;
      }
    }

    for(int i = 0; i < round; i++){
      for(int j = 0; j < num_g; j++){
        double effects = 0.0;
        for(int k = 0; k < num_g; k++){
          effects += interactions.at(j).at(k) * p_vec.at(k);
        }

        tmp_p_vec.at(j) = p_vec.at(j) + (0.5 + 0.5 * std::tanh(effects)) -
          alpha * p_vec.at(j);
      }
      p_vec = tmp_p_vec;

      for(int j = 0; j < num_g; j++){
        if(p_vec.at(j) > mid_vals.at(j) && ample.at(j) == 1 && sign.at(j) == 0){
          fluct = 1;
        }else if(p_vec.at(j) < mid_vals.at(j) && ample.at(j) == 1 && sign.at(j) == 1){
          fluct = 0;
        }
      }
    }
  }

  return(fluct);
}
