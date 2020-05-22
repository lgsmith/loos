#include "DFTMagicCircle.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;
double f(double t, vector<double> frqs) {
  double accum = 0;
  for (const auto w : frqs) {
    accum += cos(w * 2 * PI * t);
  }
  return accum;
}

int main() {
  vector<double> frqs{440, 1800, 4000};
  vector<double> signal_frqs{440, 125, 10, 4000, 320, 0};
  const double sample_frq = 8000;
  const double N_samples = 1000;
  MatrixXd sample(1, 1);
  sample << 0;
  cout << "\nSAMPLE first:\n" << sample << endl;

  DFTMagicCircle dft(sample, frqs, sample_frq, N_samples);

  for (int i = 0; i < N_samples; i++) {
    sample << f(i / sample_frq, signal_frqs);
    cout << "\nSample after " << i << "\n" << sample << "\n";
    dft(sample);
    // cout << "y1:\n";
    // for (const auto i : dft.y1)
    //   cout << i << " ";
    // cout "\ny2:\n";
    // for (const auto i : dft.y2)
    //   cout << i << " ";
    // cout << endl;

  }
  vector<MatrixXd> J = dft.spectral_density();
  for (auto i = 0; i < frqs.size(); i++)
    cout << "For frequency: " << frqs[i] << " density: " << J[i] << "\n";
  return 0;
}