#ifndef OBSERVABLE_HPP
#define OBSERVABLE_HPP

#include <cmath> // for std::sqrt
#include <string>

class observable {
public:
  observable(std::string const& name) : name_(name), count_(0), sum_(0), esq_(0) {}
  void operator<<(double x) { sum_ += x; esq_ += x * x; ++count_; }
  double mean() const { return (count_ > 0) ? (sum_ / count_) : 0.; }
  double error() const {
    return (count_ > 1) ?
      std::sqrt((esq_ / count_ - mean() * mean()) / (count_ - 1)) : 0.;
  }
  std::ostream& write(std::ostream& os) const {
    os << name_ << " = " << this->mean() << " +- " << this->error();
    return os;
  }
private:
  std::string name_;
  unsigned int count_;
  double sum_, esq_;
};

std::ostream& operator<<(std::ostream& os, observable const& obs) {
  return obs.write(os);
}
  
#endif // OBSERVABLE_HPP
