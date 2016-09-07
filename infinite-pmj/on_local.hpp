#ifndef ON_LOCAL_HPP
#define ON_LOCAL_HPP

// O(N) local flip algorithm

class on_local {
public:
  on_local(double coupling, std::vector<std::vector<int> > const& sign, double epsilon) :
    num_sites_(sign.size()), coupling_(coupling), sign_(sign), epsilon_(epsilon) {
    if (num_sites_ <= 0) {
      std::cerr << "Error: number of sites\n"; std::exit(127);
    }
    if (coupling <= 0.0) {
      std::cerr << "Error: sign of coupling\n"; std::exit(127);
    }
    for (int i = 0; i < num_sites_; ++i) {
      if (sign_[i].size() != num_sites_) {
        std::cerr << "Error: size of sign\n"; std::exit(127);
      }
    }
    for (int i = 0; i < num_sites_; ++i) {
      if (sign_[i][i] != 0) {
        std::cerr << "Error: diagonal element of sign\n"; std::exit(127);
      }
      for (int j = (i + 1); j < num_sites_; ++j) {
        if ((sign_[i][j] != 1) && (sign_[i][j] != -1)) {
          std::cerr << "Error: offdiagonal element of sign\n"; std::exit(127);
        }
        if (sign_[i][j] != sign_[j][i]) {
          std::cerr << "Error: symmetry of sign\n"; std::exit(127);
        }
      }
    }
    if (epsilon_ < 0) {
      std::cerr << "Error: sign of epsilon\n"; std::exit(127);
    }
  }
  template<typename UNIFORM_01>
  void update(double temperature, std::vector<int>& spins, UNIFORM_01& uniform_01) {
    if (temperature <= 0.0) {
      std::cerr << "Error: sign of temperature\n"; std::exit(127);
    }
    if (spins.size() != num_sites_) {
      std::cerr << "Error: size of spins\n"; std::exit(127);
    }
    double lambda = (2.0 + epsilon_) * num_sites_ * coupling_ / temperature;
    for (int s0 = 0; s0 < num_sites_; ++s0) {
      double r = 1;
      for (double t = exp_dist(uniform_01); t < lambda; t += exp_dist(uniform_01)) {
        int s1 = num_sites_ * uniform_01();
        if (s1 != s0) {
          int sigma = spins[s0] * spins[s1];
          r *= (1 + epsilon_ - sigma * sign_[s0][s1]) / (1 + epsilon_ + sigma * sign_[s0][s1]);
        }
      }
      if (uniform_01() < r / (1 + r)) {
        spins[s0] = -spins[s0];
      }
    }
  }
protected:
  template<typename UNIFORM_01>
  static double exp_dist(UNIFORM_01& uniform_01) {
    return -std::log(1 - uniform_01());
  }
private:
  int num_sites_;
  double coupling_;
  std::vector<std::vector<int> > sign_;
  double epsilon_;
};

#endif // ON_LOCAL_HPP
