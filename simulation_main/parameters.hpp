#ifndef PARAMETERS
#define PARAMETERS

class Parameters{
public:
  int pop_size;

  int num_g;
  double alpha;

  double gamma_g;
  double mu_g;
  double gamma_b;
  double mu_b;
  double gamma_c;
  double mu_c;

  double max_g;

  double sig;
  double criteria;

  int dev_period;
  int sep;
  int sensor_period;

  int num_trial;
};

#endif
