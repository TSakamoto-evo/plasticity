#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

int main(){
  int time = 100000;
  int sep = 300;
  double criteria = 1e-4;
  double alpha = 0.2;

  std::vector<double> g;
  std::vector<std::vector<double>> b;
  std::vector<double> c;

  {
    std::ifstream ifs("simulation/ind_regi.txt");
    if(!ifs){
      std::cerr << "Fail to open the file!" << std::endl;
      std::exit(1);
    }

    std::vector<std::string> lines;
    std::string line;
    
    while (getline(ifs, line)){
      std::istringstream iss(line);
      std::string tmp_list;
      std::vector<std::string> list;

      while(getline(iss, tmp_list, '\t')){
        list.push_back(tmp_list);
      }

      if(time == std::stoi(list[0])){
        int num_g = std::stoi(list[1]);
        int index = 2;

        for(int i = 0; i < num_g; i++){
          g.push_back(std::stod(list[index]));
          index++;
        }

        for(int i = 0; i < num_g; i++){
          std::vector<double> tmp;

          for(int j = 0; j < num_g; j++){
            tmp.push_back(std::stod(list[index]));
            index++;
          }

          b.push_back(tmp);
        }

        for(int i = 0; i < num_g; i++){
          c.push_back(std::stod(list[index]));
          index++;
        }
      }
    }
  }

  {
    std::cout << "g";
    for(const auto &i: g){
      std::cout << "\t" << i;
    }
    std::cout << std::endl;

    for(const auto &i1: b){
      std::cout << "b";
      for(const auto &i2: i1){
        std::cout << "\t" << i2;
      }
      std::cout << std::endl;
    }

    std::cout << "c";
    for(const auto &i: c){
      std::cout << "\t" << i;
    }
    std::cout << std::endl;
  }

  std::vector<std::vector<double>> equ;
  std::ofstream ofs("basin.txt");

  for(int i = 0; i <= sep; i++){
    for(int j = 0; j <= sep; j++){
      double x1 = 5.0 * i / sep;
      double x2 = 5.0 * j / sep;

      double convergence = 1.0;
      int gen = 0;

      while(convergence > 1e-2 * criteria || gen < 10000){
        double px1 = (1.0 - alpha) * x1 + 0.5 * 
          (1.0 + std::tanh(b.at(0).at(0) * x1 + b.at(0).at(1) * x2));
        double px2 = (1.0 - alpha) * x2 + 0.5 * 
          (1.0 + std::tanh(b.at(1).at(0) * x1 + b.at(1).at(1) * x2));
        
        convergence = std::abs(x1 - px1) + std::abs(x2 - px2);
        x1 = px1;
        x2 = px2;

        gen++;
      }

      if(equ.size() == 0){
        std::vector<double> tmp = {x1, x2};
        equ.push_back(tmp);
      }else{
        bool already_found = 0;
        for(int k = 0; k < static_cast<int>(equ.size()); k++){
          double dist = std::abs(equ.at(k).at(0) - x1) + std::abs(equ.at(k).at(1) - x2);

          if(dist < criteria){
            already_found = 1;
            ofs << 5.0 * i / sep << "\t" << 5.0 * j / sep << "\t" << k << std::endl;
            break;
          }
        }

        if(already_found == 0){
          ofs << 5.0 * i / sep << "\t" << 5.0 * j / sep << "\t" << equ.size() << std::endl;
          
          std::vector<double> tmp = {x1, x2};
          equ.push_back(tmp);
        }
      }
    }
  }

  std::ofstream ofs2("equilibria.txt");
  for(int i = 0; i < static_cast<int>(equ.size()); i++){
    ofs2 << i << "\t" << equ.at(i).at(0) << "\t" << equ.at(i).at(1) << std::endl;
  }

  std::cout << "# of equilibria: " << equ.size() << std::endl;
}
