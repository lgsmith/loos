/* internuclear-vector-corr-projections
  Writes the z projection of all selected internuclear vectors out as a
  timeseries, Optionally finding each autocorrelation as well.
*/

/*
  This file is part of LOOS.

  LOOS (Lightweight Object-Oriented Structure library)
  Copyright (c) 2020, Louis
  Department of Biochemistry and Biophysics
  School of Medicine & Dentistry, University of Rochester

  This package (LOOS) is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation under version 3 of the License.

  This package is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "loos.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#define EIGEN_USE_THREADS

using namespace std;
using namespace loos;
using namespace Eigen;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

const string fullHelpMsg = "XXX";

// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  string approxes = "types:\nISA\nSD\n";
  // clang-format off
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("time-series,t", po::value<string>(&ts)->default_value(""),
       "(Over)write time series to provided file name. If empty, then not written.")
      ("nthreads,n", po::value<int>(&t)->default_value(1),
       "Number of threads to use for the calculation.")
      ("gamma,g", po::value<double>(&gamma)->default_value(42.58),
      "Set the gyromagnetic ratio, in MHz/T.")
      ("field-strength,B", po::value<double>(&B)->default_value(14.1),
       "Set the value of the experimental field strength, in T.")
      ("sampling-freq,f", po::value<double>(&f)->default_value(1),
       "time-spacing of samples from trajectory, in GHz (frames/ns). A frequency of zero is an error.")
      ("mix,m", po::value<double>(&m)->default_value(0),
      "NOE mixing time, in milliseconds.")
      ("buildup-curve-range", po::value<std::string>(&buildup_range), 
       "Which mixing times to write out for plotting (matlab style range, overrides -m).")
      ("isa,I", po::bool_switch(&isa),
       "If thrown, report relative NOEs without any spin-relaxation.")
      ("initial-magnetization,M", po::value<double>(&M)->default_value(1),
       "Initial magnetization (M_0) at t=0. If 1 is used, All NOEs relative.");
      
  }
  // clang-format on
  // The print() function returns a string that describes what all the
  // options are set to (for logging purposes)
  string print() const {
    ostringstream oss;
    oss << boost::format(
               "ts=%s,w=%d,B=%d,f=%d,t=%d,m=%d,M=%d,buildup-range=%s") %
               ts % gamma % B % f % t % m % M % buildup_range;
    return (oss.str());
  }
  bool isa;
  string ts;
  string buildup_range;
  double gamma, B, f, m, M;
  int t;
};
// @endcond

// Class for 'magic circle oscillator' (a-la Clay Turner) selective DFTs
template <typename SampleType, typename FrqType> class DFTMagicCircle {
private:
  std::vector<Eigen::MatrixBase<SampleType>> y1;
  std::vector<Eigen::MatrixBase<SampleType>> y2;
  std::vector<FrqType> K;
  std::vector<Eigen::MatrixBase<SampleType>> J;

public:
  // needs frequencies in Hertz
  DFTMagicCircle(Eigen::MatrixBase<SampleType> &empty_sample,
                 std::vector<FrqType> &frqs, double sampling_frq,
                 double n_samples);
  void operator()(const Eigen::MatrixBase<SampleType> &sample);
  std::vector<Eigen::MatrixBase<SampleType>> spectral_density(void);
  ~DFTMagicCircle();
};

template <typename SampleType, typename FrqType>
DFTMagicCircle<SampleType, FrqType>::DFTMagicCircle(
    Eigen::MatrixBase<SampleType> &empty_sample, std::vector<FrqType> &frqs,
    double sampling_frq, double n_samples) {
  // buildup k vector with sinusoids corresponding to tracked frqs.
  for (const auto f : frqs) {
    // convert to radians per sample over two, take sin, then multiply by 2
    K.push_back(2 * std::sin((PI / n_samples) * std::floor(f / sampling_frq)));
    // for each frq, set up both recursion half-step states to zero.
    y1.push_back(Eigen::MatrixBase<SampleType>::Zero(empty_sample.rows(),
                                                     empty_sample.cols()));
    y2.push_back(Eigen::MatrixBase<SampleType>::Zero(empty_sample.rows(),
                                                     empty_sample.cols()));
  }
}

template <typename SampleType, typename FrqType>
inline void DFTMagicCircle<SampleType, FrqType>::
operator()(const Eigen::MatrixBase<SampleType> &sample) {
  for (auto i = 0; i < K.size(); i++) {
    // do both DFT sinusoid half steps now
    y2[i].template triangularView<Eigen::Lower>() += sample - K[i] * y1[i];
    y1[i].template triangularView<Eigen::Lower>() += K[i] * y2[i];
  }
}

template <typename SampleType, typename FrqType>
inline std::vector<Eigen::MatrixBase<SampleType>>
DFTMagicCircle<SampleType, FrqType>::spectral_density() {
  for (auto i = 0; i < K.size(); i++) {
    J.push_back(
      y1[i].cwiseAbs2() + y2[i].cwiseAbs2() - (K[i] * y1[i].cwiseProduct(y2[i]))
    );
  }
  return J;
}

template <typename SampleType, typename FrqType>
DFTMagicCircle<SampleType, FrqType>::~DFTMagicCircle() {}

// time conversions
const double ghz2Hz = 1e9;
const double mhz2Hz = 1e6;
const double ms2s = 1e-3;

// second order spherical harmonic, with multiples factored out.
inline const double Y_2_0(GCoord a, GCoord b) {
  const GCoord diff = a - b;
  const double length = diff.length();
  const double r3 = length * length * length;
  return 1.5 * diff[2] * diff[2] / (r3 * length * length) - (0.5 / r3);
}

int main(int argc, char *argv[]) {

  string header = invocationHeader(argc, argv);

  opts::BasicOptions *bopts = new opts::BasicOptions;
  opts::BasicSelection *sopts = new opts::BasicSelection;
  opts::MultiTrajOptions *mtopts = new opts::MultiTrajOptions;
  ToolOptions *topts = new ToolOptions;

  // combine options
  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mtopts).add(topts);

  // Parse the command-line.  If an error occurred, help will already
  // be displayed and it will return a FALSE value.
  if (!options.parse(argc, argv))
    exit(-1);

  // Pull the model from the options object (it will include coordinates)
  AtomicGroup model = mtopts->model;

  // Pull out the trajectory...
  pTraj traj = mtopts->trajectory;
  if (mtopts->mtraj.nframes() == 0) {
    cout << "ERROR: can't work with zero frames.\n";
    exit(-1);
  }

  // Select the desired atoms to operate over...
  AtomicGroup nuclei = selectAtoms(model, sopts->selection);
  const Eigen::Index N = (Eigen::Index)nuclei.size();
  vector<vector<uint>> refindex = {{2, 3}};
  const double refdist = nuclei[refindex[0][0]]->coords().distance(
      nuclei[refindex[0][1]]->coords());

  // NMR precalculations
  const double mu0 =
      1.25663706212e4; // wikipedia, H/Angstrom (10^10 * value in H/m)
  const double gamma =
      topts->gamma * 2 * PI * mhz2Hz;  // convert Gamma from mHz/T to Rad/s*T
  const double hbar = 1.054571817e-34; // wikipedia, J*s
  const double N_A = 6.02214076e24;    // Wikipedia, Avogadro's Constant
  // dipolar interaction constant, unit distance per Mole
  const double dd = gamma * gamma * mu0 * hbar / (4 * PI);
  const double dd2 = dd * dd * N_A * 5.0 / (PI * 16);

  // Magic circle oscillator precomputation:
  // Larmor frequency, in Hz
  const double omega = gamma * topts->B * mhz2Hz;
  // convert to bin number omega
  const double omega_bin = floor(omega / (topts->f * ghz2Hz));
  // for magic circle, the customary 2 Pi bin_num/nSamples is divided by 2 to
  // produce two 90 degree half-step updates per sample from the time-series.
  const double omega_rad_sample = omega_bin * PI / mtopts->frameList().size();
  // DFT bin corresponding to the above omega
  // we need three frequencies; 0, omega, and two times omega
  const double k = 2 * std::sin(omega_rad_sample);
  const double k2 = 2 * std::sin(2 * omega_rad_sample);
  // matricies to hold the intermediate values we need for the three frequencies
  // MatrixXd y1zero = MatrixXd::Zero(N, N); // unneeded, always zero for k = 0.
  MatrixXd y2zero = MatrixXd::Zero(N, N);
  MatrixXd y1omega = MatrixXd::Zero(N, N);
  MatrixXd y2omega = MatrixXd::Zero(N, N);
  MatrixXd y1omega2 = MatrixXd::Zero(N, N);
  MatrixXd y2omega2 = MatrixXd::Zero(N, N);

  MatrixXd y_2_0 = MatrixXd::Zero(N, N);
  GCoord diff(0, 0, 0);
  // Now iterate over all frames in the skipped & strided trajectory
  for (auto f : mtopts->frameList()) {

    // Update the coordinates ONLY for the subset
    traj->readFrame(f);
    traj->updateGroupCoords(nuclei);
    // Get coords into a tensor (presuming unrolling below)
    for (auto i = 0; i < N; i++) {
      for (auto j = 0; j < i; j++) {
        y_2_0(i, j) = Y_2_0(nuclei[i]->coords(), nuclei[j]->coords());
      }
    }
    // cout << "this is frame: " << f << "\n"
    //  << y_2_0 << "\n";
    // cout << cos_r3 << "\n\n";
    // first update the y2s for the two frequencies
    // since k-value is zero for zero freq, this is just accumulating cosine
    // term.
    y2zero.triangularView<Lower>() += y_2_0;
    // these two use k and k2 defined above
    y2omega.triangularView<Lower>() += y_2_0 - k * y1omega;

    y2omega2.triangularView<Lower>() += y_2_0 - k2 * y1omega2;
    // now compute the y1 from these y2, noting that y1zero will be zero
    y1omega.triangularView<Lower>() += k * y2omega;
    y1omega2.triangularView<Lower>() += k2 * y2omega2;
  }
  // compute the spectral densities at each of the three needed freqs
  // following formula E = y1**2 +y2**2 - k*y1*y2
  // Note that here all calculations should be array calx, not matrix calx.
  MatrixXd Jzero = y2zero.cwiseAbs2();
  cout << "This is Jzero:\n" << Jzero << "\n";
  MatrixXd Jomega = y1omega.cwiseAbs2() + y2omega.cwiseAbs2() -
                    k * y1omega.cwiseProduct(y2omega);

  cout << "This is Jomega:\n" << Jomega << "\n";

  MatrixXd Jomega2 = y1omega2.cwiseAbs2() + y2omega2.cwiseAbs2() -
                     k2 * y1omega2.cwiseProduct(y2omega2);
  cout << "This is Jomega2:\n" << Jomega2 << "\n";
  // Comput sigma_{ij} and rho_i following Chalmers et al.
  // Sigma is the cross-relaxation rate, and is
  // the sum over the full power and omega spectral densities.
  // MatrixXd sigma = 6 * Jomega2 - Jzero;
  MatrixXd rho((Jzero + (3 * Jomega) + (6 * Jomega2)).selfadjointView<Lower>());
  MatrixXd R(6 * Jomega2 - Jzero);
  R.diagonal(0) = rho.colwise().sum();
  cout << "this is R with sigma added, but not in correct units:\n";
  cout << R << endl;
  R *= dd2;
  cout << "this is R:\n";
  cout << R << endl;
  // if (topts->isa) {
  // do report based on ISA

  // } else {
  SelfAdjointEigenSolver<MatrixXd> es(R);
  cout << es.eigenvalues() << endl;
  MatrixXd evolved_vals = (es.eigenvalues() * (-topts->m * ms2s))
                              .array()
                              .exp()
                              .matrix()
                              .asDiagonal();
  cout << "this is evolved_vals:\n" << evolved_vals << endl;
  MatrixXd intensities = es.eigenvectors() * evolved_vals *
                         es.eigenvectors().inverse() *
                         (topts->M * MatrixXd::Identity(N, N));
  cout << intensities << endl;
  // }
  // create tab delimited intensity report, below
  cout << "reference intensity and distance:\n"
       << intensities(refindex[0][0], refindex[0][1]) << " " << refdist << "\n";
  cout << "# resname\tresid\tname\tindex\tresname\tresid\tname\tindex\tvol\tme"
          "an_r";
  for (auto i = 0; i < N; i++) {
    pAtom ith = nuclei[i];
    for (auto j = i + 1; j < N; j++) {
      pAtom jth = nuclei[j];
      cout << "\n"
           << ith->resname() << "\t" << ith->resid() << "\t" << ith->name()
           << "\t" << ith->index() << "\t" << jth->resname() << "\t"
           << jth->resid() << "\t" << jth->name() << "\t" << jth->index()
           << "\t" << intensities(i, j) << "\t"
           << pow(intensities(i, j) /
                      intensities(refindex[0][0], refindex[0][1]),
                  -1.0 / 6) *
                  refdist;
    }
  }

  ComputationInfo es_info = es.info();
  if (es_info == Success)
    cout << "\nEigendecomposition successful.\n";
  if (es_info == NumericalIssue)
    cout << "\nEigendecomposition ran into a numerical issue.\n";
  if (es_info == NoConvergence)
    cout << "\nEigendecomposition did not converge.\n";
  if (es_info == InvalidInput)
    cout << "\nEigendecomposition was given invalid input.\n";
}
