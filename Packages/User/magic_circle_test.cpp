#include "DFTMagicCircle.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;
double f(double t, vector<double> &frqs) {
  double accum = 0;
  for (const auto w : frqs) {
    accum += cos(w * 2 * M_1_PI * t);
  }
  return accum;
}

int main() {
  vector<double> frqs{2000, 1800, 4000};
  vector<double> signal_frqs{2000, 125, 10, 4000, 320, 0};
  const double sample_frq = 8000;
  const int N_samples = 16000;
  MatrixXd sample(2, 2);
  sample << 0, 0, 0, 0;
  vector<MatrixXd> samples;
  cout << "\nSAMPLE first:\n" << sample;

  for (int i = 0; i < N_samples; i++) {
    sample(1, 0) = f(i / sample_frq, signal_frqs);
    samples.push_back(sample);
  }

  DFTMagicCircle dft(sample, frqs, sample_frq, N_samples);
  for (auto s : samples) {
    cout << "\n" << s << "\n";
    dft(s);
    for (uint i = 0; i < frqs.size(); i++) {
      cout << dft.get_y1()[i] << " " << dft.get_y2()[i] << " ";
    }
  }
  vector<MatrixXd> J = dft.spectral_density();
  // gain for DFT is approx N/2, so gain on spectrum is (N/2)^2
  const double gain = N_samples * N_samples / 4;
  for (uint i = 0; i < frqs.size(); i++) {
    cout << "\nFor frequency:\n"
         << frqs[i] << "  " <<  dft.get_k()[i] << "\ndensity:\n"
         << J[i] / gain << "\n";
  }
  return 0;
}