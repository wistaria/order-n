// Copyright (C) 2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

#ifndef NAIVE_MC_HPP
#define NAIVE_MC_HPP

class naive_mc {
public:
  naive_mc(double coupling, std::vector<std::vector<int> > const& sign) :
    num_sites_(sign.size()), coupling_(coupling), sign_(sign) {
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
  }
  template<typename UNIFORM_01>
  void update(double temperature, std::vector<int>& spins, UNIFORM_01& uniform_01) {
    if (temperature <= 0.0) {
      std::cerr << "Error: sign of temperature\n"; std::exit(127);
    }
    if (spins.size() != num_sites_) {
      std::cerr << "Error: size of spins\n"; std::exit(127);
    }
    double beta = 1.0 / temperature;
    for (int i = 0; i < num_sites_; ++i) {
      double diff = 0;
      for (int j = 0; j < num_sites_; ++j) {
        diff += sign_[i][j] * spins[i] * spins[j];
      }
      if (uniform_01() < 0.5 * (1 + std::tanh(-beta * coupling_ * diff))) {
        spins[i] = -spins[i];
      }
    }
  }
private:
  int num_sites_;
  double coupling_;
  std::vector<std::vector<int> > sign_;
};

#endif // NAIVE_MC_HPP
