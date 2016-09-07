#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/timer.hpp>
#include <observable.hpp>
#include <power.hpp>
#include "on_local.hpp"

int main(int argc, char** argv) {
  unsigned int seed;
  unsigned int N;
  double prob; // probability of interaction is ferromagnetic
  double temperature;
  unsigned int sweeps, thermalization;
  double epsilon;
  if (argc >=7) {
    seed = boost::lexical_cast<int>(argv[1]);
    N = boost::lexical_cast<double>(argv[2]);
    prob = boost::lexical_cast<double>(argv[3]);
    temperature = boost::lexical_cast<double>(argv[4]);
    sweeps = boost::lexical_cast<int>(argv[5]);
    epsilon = boost::lexical_cast<double>(argv[6]);
  } else {
    std::cin >> seed >> N >> prob >> temperature >> sweeps >> epsilon;
  }
  thermalization = (sweeps >> 3);
  std::cout << "Order-N MC for infinite-range +/- J model" << std::endl
            << "seed           = " << seed << std::endl
            << "N              = " << N << std::endl
            << "P_ferro        = " << prob << std::endl
            << "temperature    = " << temperature << std::endl
            << "sweeps         = " << sweeps << std::endl
            << "thermalization = " << thermalization << std::endl
            << "epsilon        = " << epsilon << std::endl;

  // random number generators
  boost::mt19937 eng(seed);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    uniform_01(eng, boost::uniform_real<>());

  // coupling constant
  double coupling = 1.0 / N;
  std::vector<std::vector<int> > sign(N, std::vector<int>(N, 0));
  for (int i = 0; i < N; ++i) {
    for (int j = i + 1; j < N; ++j) {
      sign[i][j] = sign[j][i] = (uniform_01() < prob ? 1 : -1);
    }
  }

  // configuration
  std::vector<int> spins(N, 1); // all up
  
  // MC updater
  on_local mc(coupling, sign, epsilon);

  // observables
  observable magnetization2("Magnetization^2");
  observable magnetization4("Magnetization^4");
  
  // MC step
  boost::timer tm;
  for (unsigned int mcs = 0; mcs < thermalization + sweeps; ++mcs) {
    mc.update(temperature, spins, uniform_01);
    if (mcs >= thermalization) {
      double mag = 0;
      for (int i = 0; i < N; ++i) mag += spins[i];
      magnetization2 << power2(mag);
      magnetization4 << power4(mag);
    }
  }

  double elapsed = tm.elapsed();
  std::clog << "Elapsed time = " << elapsed << " sec\n"
            << "Speed = " << (thermalization + sweeps) / elapsed << " MCS/sec\n";
  std::cout << magnetization2 << std::endl
            << magnetization4 << std::endl;
}
