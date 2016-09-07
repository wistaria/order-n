// Copyright (C) 2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

#include <boost/lexical_cast.hpp>

#ifndef OPTIONS_HPP
#define OPTIONS_HPP

struct options {
  unsigned int seed;
  unsigned int N;
  double p_ferro; // probability of interaction is ferromagnetic
  double temperature;
  unsigned int sweeps, thermalization;
  double epsilon;

  options(int argc, char** argv) : epsilon(0.1) {
    if (argc >= 6) {
      seed = boost::lexical_cast<int>(argv[1]);
      N = boost::lexical_cast<double>(argv[2]);
      p_ferro = boost::lexical_cast<double>(argv[3]);
      temperature = boost::lexical_cast<double>(argv[4]);
      sweeps = boost::lexical_cast<int>(argv[5]);
      if (argc >= 7) epsilon = boost::lexical_cast<double>(argv[6]);
    } else {
      std::cin >> seed >> N >> p_ferro >> temperature >> sweeps >> epsilon;
    }
    thermalization = (sweeps >> 3);
    std::cout << "seed               = " << seed << std::endl
              << "N                  = " << N << std::endl
              << "P_ferro            = " << p_ferro << std::endl
              << "temperature        = " << temperature << std::endl
              << "sweeps             = " << sweeps << std::endl
              << "thermalization     = " << thermalization << std::endl
              << "epsilon (for O(N)) = " << epsilon << std::endl;
  }
};

#endif // OPTIONS_HPP
