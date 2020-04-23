/*
  internuclear-vector-corr-projections
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
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/CXX11/ThreadPool>

using namespace std;
using namespace loos;
using namespace Eigen;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;

const string fullHelpMsg = "XXX";

// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  // clang-format off
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("time-series,t", po::value<string>(&ts)->default_value(""),
       "(Over)write time series to provided file name. If empty, then not written.")
      ("nthreads,n", po::value<int>(&t)->default_value(1),
       "Number of threads to use for the calculation.")
      ("gamma,g", po::value<double>(&w)->default_value(42.58),
      "Set the larmor frequency to arg.")
      ("field-strength,B", po::value<double>(&B)->default_value(14.1),
       "Set the value of the experimental field strength.")
      ("sampling-freq,f", po::value<double>(&f)->default_value(1),
       "time-spacing of samples from trajectory, in GHz (frames/ns). A frequency of zero is an error.")
      ("mix,m", po::value<double>(&m)->default_value(0),
      "Mixing time for experiment. A value of zero means no SD.")
      ("initial-magnetization,M", po::value<double>(&M)->default_value(1),
       "Initial magnetization (M_0) at t=0. If 1 is used, All NOEs relative.")
      ;
  }
  // clang-format on
  // The print() function returns a string that describes what all the
  // options are set to (for logging purposes)
  string print() const {
    ostringstream oss;
    oss << boost::format("ts=%s,w=%d,B=%d,f=%d,t=%d,m=%d,M=%d") % ts % w % B %
               f % t % m % M;
    return (oss.str());
  }
  string ts;
  double w, B, f, m, M;
  int t;
};
// @endcond

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
  const auto N = nuclei.size();

  // NMR precalculations
  double dd = 1; // dipolar interaction constant, unit distance

  // Broadcast dims for multiplying the recurrence tensors by K
  Eigen::array<int, 3> bcRecurrence({3, N, N});
  // Magic circle oscillator precomputation
  double omega = topts->w * topts->B;
  // we need three frequencies; 0, omega, and two times omega
  double k = sin(PI * omega / mtopts->mtraj.nframes());
  double k2 = sin(PI * 2 * omega / mtopts->mtraj.nframes());
  TensorFixedSize<double, Sizes<3>, RowMajor> K_base;
  Tensor<double, 3, RowMajor> K(3, N, N);
  // constant sinusoids to be generated.
  K_base.setValues({0, k, k2});
  K = K_base.broadcast(bcRecurrence);

  // setup for eigen tensor ops below
  Eigen::array<int, 2> bcChip({N, 3});
  Tensor<double, 2, RowMajor> coords(N, 3);
  Tensor<double, 3, RowMajor> internuclear_vectors(N, N, 3);
  // Define tensors to hold recurrence relation values for magic circle
  // oscillator.
  Tensor<double, 3, RowMajor> p1(3, N, N);
  Tensor<double, 3, RowMajor> p2(3, N, N);
  p1.setZero();
  p2.setZero();

  // Dimension object to pick out the coordinates from the internuclear vectors
  // tensor for summing
  Eigen::array<int, 1> dimCoords({2});

  // threading setup for tensor ops
  ThreadPool pool(topts->t);
  ThreadPoolDevice threader(&pool, topts->t);
  // threading for openmp parallelization of eigen matrix/vector products
  setNbThreads(topts->t);

  // Now iterate over all frames in the skipped & strided trajectory
  while (traj->readFrame()) {

    // Update the coordinates ONLY for the subset
    traj->updateGroupCoords(nuclei);
    // Get coords into a tensor (presuming unrolling below)
    for (auto i = 0; i < nuclei.size(); i++)
      for (auto j = 0; j < 3; j++)
        coords(i, j) = nuclei[i]->coords()[j];

    // compute interatomic vectors
    for (auto i = 0; i < N; i++) {
      internuclear_vectors.chip(i, 0).device(threader) =
          coords - coords.chip(i, 0).eval().broadcast(bcChip);
    }
    cout << "\ninvs from frame : " << mtopts->mtraj.currentFrame() << "\n"
         << internuclear_vectors << endl;
    // get z coords, divided by length of inv^4 (cosine(theta)/r^3); auto causes
    // lazy eval.
    auto Xs = internuclear_vectors.chip(2, 2) /
              (internuclear_vectors.square().sum(dimCoords).sqrt()).pow(4);
    // compute magic circle oscillator recurrence relations for this frame
    p2.device(threader) = p2 - K * p1 + Xs.eval().broadcast(bcRecurrence);
    p1.device(threader) = p1 + K * p2;
  }
  // this expression records the squared value of the spectral density at the
  // three freqs.
  auto J_base = p1 * p1 + p2 * p2 - K * p1 * p2;

  Tensor<double, 3, RowMajor> nans(3, N, N);
  double not_a_number = std::numeric_limits<double>::quiet_NaN();
  nans.setConstant(not_a_number);
  cout << nans << endl;
  Tensor<double, 3, RowMajor> zeros(3, N, N);
  zeros.setZero();
  auto J = (J_base != J_base).select(zeros, J_base);
Tensor<double, 3, RowMajor> J_t = J;
  cout << "J tensor:\n" << J_t << endl;
  // Comput sigma_{ij} and rho_i following Chalmers et al.
  // Sigma is the cross-relaxation rate, and is
  // the sum over the full power and omega spectral densities.
  auto sigma = J.chip(0, 0) - 6 * J.chip(2, 0);
  // rho is the diagonal of the relaxation matrix,
  // and is the sum over all non-diagonal elements.
  auto rho = (J.chip(0, 0) + 3 * J.chip(1, 0) + 6 * J.chip(2, 0))
                 .sum(Eigen::array<int, 1>({1}));

  Tensor<double, 2, RowMajor> R_t(N, N);
  R_t.device(threader) = sigma + rho;
  cout << "R tensor:\n" << R_t << endl;
  auto R = MatrixXd::Map(R_t.data(), N, N);
  SelfAdjointEigenSolver<MatrixXd> es(R);
  MatrixXd evolved_evs =
      (es.eigenvalues() * topts->m).array().exp().matrix().asDiagonal();
  MatrixXd intensities = es.eigenvectors() * evolved_evs *
                         es.eigenvectors().inverse() *
                         (topts->M * MatrixXd::Identity(N, N));
  cout << intensities << endl;
}
