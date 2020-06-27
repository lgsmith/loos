#include <Eigen/Dense>
#include <cmath>
#include <vector>

typedef Eigen::MatrixXd SampleType;
class DFTMagicCircle {
private:
  std::vector<SampleType> J;
  std::vector<SampleType> y2; // needs frequencies in Hertz
  std::vector<SampleType> y1;
  std::vector<double> k;
  double scale;

public:
  // needs frequencies in Hertz
  DFTMagicCircle(SampleType &empty_sample, const std::vector<double> &frqs,
                 const double fs, const unsigned int n_samples);
  void operator()(const SampleType &sample);
  std::vector<SampleType> power_spectral_density(bool recompute = false);
  // setters and getters
  std::vector<SampleType> get_y2(void);
  void set_y2(std::vector<SampleType> &recur2);
  std::vector<SampleType> get_y1(void);
  void set_y1(std::vector<SampleType> &recur1);
  std::vector<double> get_k(void);
  void set_k(std::vector<double> &frqs);
  // end setters and getters
  ~DFTMagicCircle();
};

DFTMagicCircle::DFTMagicCircle(SampleType &empty_sample,
                               const std::vector<double> &frqs,
                               const double fs,
                               const unsigned int n_samples) {
  // compute scale-factor now, for use returning power spectral density
  scale = 1.0 / (static_cast<double>(n_samples) * fs);
  // buildup k vector with sinusoids corresponding to tracked frqs.
  for (const auto f : frqs) {
    // convert to radians per sample over two, take sin, then multiply by 2
    k.push_back(2 * std::sin(f * M_1_PI / fs)); // note pi, not 2*pi
    // for each frq, set up both recursion half-step states to zero.
    y1.push_back(SampleType::Zero(empty_sample.rows(), empty_sample.cols()));
    y2.push_back(SampleType::Zero(empty_sample.rows(), empty_sample.cols()));
  }
}
// begin getters and setters
// recurrence element 2
std::vector<SampleType> DFTMagicCircle::get_y2(void) { return y2; }
void DFTMagicCircle::set_y2(std::vector<SampleType> &recur2) { y2 = recur2; }
// recurrence element 1
std::vector<SampleType> DFTMagicCircle::get_y1(void) { return y1; }
void DFTMagicCircle::set_y1(std::vector<SampleType> &recur1) { y1 = recur1; }
// sinusoids projected onto
std::vector<double> DFTMagicCircle::get_k(void) { return k; }
void DFTMagicCircle::set_k(std::vector<double> &frqs) { k = frqs; }
// end getters and setters
inline void DFTMagicCircle::operator()(const SampleType &sample) {
  for (uint i = 0; i < k.size(); i++) {
    // do both DFT sinusoid half steps now
    y2[i].triangularView<Eigen::Lower>() += sample - k[i] * y1[i];
    y1[i].triangularView<Eigen::Lower>() += k[i] * y2[i];
  }
}
// power spectral density
inline std::vector<SampleType>
DFTMagicCircle::power_spectral_density(bool recompute) {
  if (J.size() != k.size() || recompute) {
    for (uint i = 0; i < k.size(); i++) {
      J.push_back(
        scale * (y1[i].cwiseAbs2() + y2[i].cwiseAbs2() -
                  (k[i] * y1[i].cwiseProduct(y2[i])))
      );
    }
  }

  return J;
}

DFTMagicCircle::~DFTMagicCircle() {}
