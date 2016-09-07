// Copyright (C) 2016 by Synge Todo <wistaria@phys.s.u-tokyo.ac.jp>

#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <observable.hpp>
#include <power.hpp>
#include "options.hpp"
#include "order_n.hpp"

int main(int argc, char** argv) {
  std::cout << "Order-N MC for infinite-range +/- J model" << std::endl;
  options opt(argc, argv);

  // random number generators
  boost::mt19937 eng(opt.seed);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    uniform_01(eng, boost::uniform_real<>());

  // coupling constant
  double coupling = 1.0 / opt.N;
  std::vector<std::vector<int> > sign(opt.N, std::vector<int>(opt.N, 0));
  for (int i = 0; i < opt.N; ++i) {
    for (int j = i + 1; j < opt.N; ++j) {
      sign[i][j] = sign[j][i] = (uniform_01() < opt.p_ferro ? 1 : -1);
    }
  }

  // configuration
  std::vector<int> spins(opt.N, 1); // all up
  
  // MC updater
  order_n mc(coupling, sign, opt.epsilon);

  // observables
  observable magnetization2("Magnetization^2");
  observable magnetization4("Magnetization^4");
  
  // MC step
  boost::timer tm;
  for (unsigned int mcs = 0; mcs < opt.thermalization + opt.sweeps; ++mcs) {
    mc.update(opt.temperature, spins, uniform_01);
    if (mcs >= opt.thermalization) {
      double mag = 0;
      for (int i = 0; i < opt.N; ++i) mag += spins[i];
      magnetization2 << power2(mag);
      magnetization4 << power4(mag);
    }
  }

  double elapsed = tm.elapsed();
  std::clog << "Elapsed time = " << elapsed << " sec\n"
            << "Speed = " << (opt.thermalization + opt.sweeps) / elapsed << " MCS/sec\n";
  std::cout << magnetization2 << std::endl
            << magnetization4 << std::endl;
}
