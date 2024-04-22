#ifndef PARAMETERS
#define PARAMETERS

class Parameters{
public:
  int pop_size;

  int num_g;
  double alpha;

  double mu_g;
  double p_g;
  double mu_b;
  double p_b;
  double mu_c;
  double p_c;

  double max_g;

  double sig;
  double criteria;

  int dev_period;
  int sep;
};

#endif
