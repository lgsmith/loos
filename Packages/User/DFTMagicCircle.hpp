#include <Eigen/Dense>
#include <cmath>
#include <vector>

const double PI = M_PI;

class DFTMagicCircle {
private:
  std::vector<Eigen::MatrixXd> J;

public:
  // needs frequencies in Hertz
  DFTMagicCircle(Eigen::MatrixXd &empty_sample, const std::vector<double> &frqs,
                 const double sampling_rate, const long unsigned int n_samples);
  void operator()(const Eigen::MatrixXd &sample);
  std::vector<Eigen::MatrixXd> y2; // needs frequencies in Hertz
  std::vector<Eigen::MatrixXd> y1;
  std::vector<double> K;
  std::vector<Eigen::MatrixXd> spectral_density(void);
  ~DFTMagicCircle();
};

DFTMagicCircle::DFTMagicCircle(Eigen::MatrixXd &empty_sample,
                               const std::vector<double> &frqs,
                               const double sampling_rate,
                               const long unsigned int n_samples) {
  // buildup k vector with sinusoids corresponding to tracked frqs.
  for (const auto f : frqs) {
    // convert to radians per sample over two, take sin, then multiply by 2
    K.push_back(2 * std::sin((PI / n_samples) * std::floor(f / sampling_rate)));
    // for each frq, set up both recursion half-step states to zero.
    y1.push_back(
        Eigen::MatrixXd::Zero(empty_sample.rows(), empty_sample.cols()));
    y2.push_back(
        Eigen::MatrixXd::Zero(empty_sample.rows(), empty_sample.cols()));
  }
}

inline void DFTMagicCircle::operator()(const Eigen::MatrixXd &sample) {
  for (auto i = 0; i < K.size(); i++) {
    // do both DFT sinusoid half steps now
    y2[i].triangularView<Eigen::Lower>() += sample - K[i] * y1[i];
    y1[i].triangularView<Eigen::Lower>() += K[i] * y2[i];
  }
}

inline std::vector<Eigen::MatrixXd> DFTMagicCircle::spectral_density() {
  for (auto i = 0; i < K.size(); i++) {
    J.push_back(y1[i].cwiseAbs2() + y2[i].cwiseAbs2() -
                (K[i] * y1[i].cwiseProduct(y2[i])));
  }
  return J;
}

DFTMagicCircle::~DFTMagicCircle() {}
