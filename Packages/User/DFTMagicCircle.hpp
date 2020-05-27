#include <Eigen/Dense>
#include <cmath>
#include <vector>

typedef Eigen::MatrixXd SampleType;
class DFTMagicCircle {
private:
  std::vector<SampleType> J;

public:
  // needs frequencies in Hertz
  DFTMagicCircle(SampleType &empty_sample, const std::vector<double> &frqs,
                 const double sampling_rate, const unsigned int n_samples);
  void operator()(const SampleType &sample);
  std::vector<SampleType> y2; // needs frequencies in Hertz
  std::vector<SampleType> y1;
  std::vector<double> k;
  std::vector<SampleType> spectral_density(void);
  ~DFTMagicCircle();
};

DFTMagicCircle::DFTMagicCircle(SampleType &empty_sample,
                               const std::vector<double> &frqs,
                               const double sampling_rate,
                               const unsigned int n_samples) {
  // buildup k vector with sinusoids corresponding to tracked frqs.
  for (const auto f : frqs) {
    // convert to radians per sample over two, take sin, then multiply by 2
    k.push_back(2 * std::sin(f * M_1_PI / sampling_rate)); // note pi, not 2*pi
    // for each frq, set up both recursion half-step states to zero.
    y1.push_back(
        SampleType::Zero(empty_sample.rows(), empty_sample.cols()));
    y2.push_back(
        SampleType::Zero(empty_sample.rows(), empty_sample.cols()));
  }
}

inline void DFTMagicCircle::operator()(const SampleType &sample) {
  for (auto i = 0; i < k.size(); i++) {
    // do both DFT sinusoid half steps now
    y2[i].triangularView<Eigen::Lower>() += sample - k[i] * y1[i];
    y1[i].triangularView<Eigen::Lower>() += k[i] * y2[i];
  }
}

inline std::vector<SampleType> DFTMagicCircle::spectral_density() {
  for (auto i = 0; i < k.size(); i++) {
    J.push_back(y1[i].cwiseAbs2() + y2[i].cwiseAbs2() -
                (k[i] * y1[i].cwiseProduct(y2[i])));
  }
  return J;
}

DFTMagicCircle::~DFTMagicCircle() {}
