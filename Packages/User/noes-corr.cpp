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

#include "DFTMagicCircle.hpp"
#include "dft_helpers.hpp"
#include "loos.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <string>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>
#define EIGEN_USE_THREADS

using namespace std;
using namespace loos;
using namespace Eigen;

namespace opts = loos::OptionsFramework;
namespace po = loos::OptionsFramework::po;
namespace pt = boost::property_tree;

const string fullHelpMsg = "XXX";

// @cond TOOLS_INTERNAL
class ToolOptions : public opts::OptionsPackage {
public:
  // clang-format off
  void addGeneric(po::options_description& o) {
    o.add_options()
      ("gamma,g", po::value<double>(&gamma)->default_value(42.577478518),
       "Set the gyromagnetic ratio, in MHz/T.")
      ("field-strength,B", po::value<double>(&B)->default_value(14.1),
       "Set the value of the experimental field strength, in T.")
      ("sampling-freq,f", po::value<double>(&f)->default_value(1),
       "time-spacing of samples from trajectory, in GHz (frames/ns). A frequency of zero is an error.")
      ("mix,m", po::value<double>(&m)->default_value(0),
       "NOE mixing time, in milliseconds.")
      ("tau", po::value<double>(&tau)->default_value(0),
       "set a rigid-body correlation time, in ns.")
      ("buildup-curve-range", po::value<std::string>(&buildup_range), 
       "Which mixing times to write out for plotting (matlab style range, overrides -m).")
      ("bin-width", po::value<double>(&bin_width)->default_value(10.0),
       "Width for DFT bins, in kHz.")
      ("isa,I", po::bool_switch(&isa)->default_value(false),
       "If thrown, report relative NOEs without any spin-relaxation.")
      ("spectral-density,J", po::bool_switch(&spectral_density)->default_value(false),
       "If thrown, use spectral densities only at 0, w, and 2w for relaxation matrix elts.")
      ("correlation,C", po::bool_switch(corr)->default_value(false),
       "if thrown, compute full DFT of result; Output both spectral-density NOEs and correlation functions.")
      ("initial-magnetization,M", po::value<double>(&M)->default_value(1),
       "Initial magnetization (M_0) at t=0.");
      
  }
  // clang-format on
  // The print() function returns a string that describes what all the
  // options are set to (for logging purposes)
  string print() const {
    ostringstream oss;
    oss << boost::format("ts=%s,w=%d,B=%d,f=%d,t=%d,m=%d,M=%d,buildup_range=%s,"
                         "spectral_density=%b,isa=%b,corr=%b") %
               ts % gamma % B % f % threads % m % M % buildup_range %
               spectral_density % isa % corr;
    return (oss.str());
  }
  // this function post-processes conditional input for the tool options
  bool postConditions(po::variables_map &vm) {
    if (!buildup_range.empty()) {
      if (isa) {
        cerr << "Usage Error:\n"
             << "You've given a buildup range while also requesting\n"
             << "an independend spin-pair analysis, where there is no SD.\n";
        return false;
      } else {
        buildups = parseRange<double>(buildup_range);
        return true;
      }
    }
    if (spectral_density && corr) {
      cerr << "Usage Error: --spectral-density && --correlation\n"
           << "You've asked for two different methods of computing the "
              "spectral density. Pick one.\n";
    }
    return true;
  }
  bool isa, spectral_density, corr;
  string ts;
  string buildup_range;
  vector<double> buildups;
  double gamma, B, f, m, M, tau, overlap, bin_width;
  int threads;
};
// @endcond
// function to write evolution of TS as JSON
pt::ptree make_NOE_json(vector<MatrixXd> &intensities, vector<double> &times,
                        AtomicGroup &nuclei) {
  pt::ptree root; // we'll build this up, then return it.
  pt::ptree times_node;
  pt::ptree noes_node;
  for (const auto &time : times) {
    pt::ptree time_node;
    time_node.put("", time);
    times_node.push_back(make_pair("", time_node));
  }
  root.add_child("mixing times", times_node);
  for (size_t i = 0; i < nuclei.size(); i++) {
    pAtom ith = nuclei[i];
    for (size_t j = i + 1; j < nuclei.size(); j++) {
      pAtom jth = nuclei[j];
      ostringstream tag;
      tag << ith->resname() << ith->resid() << ith->name() << ":"
          << jth->resname() << jth->resid() << jth->name();
      pt::ptree noe_node;
      for (uint t = 0; t < times.size(); t++) {
        pt::ptree vol_node;
        vol_node.put("", intensities[t](i, j));
        noe_node.push_back(make_pair("", vol_node));
      }
      noes_node.add_child(tag.str(), noe_node);
    }
  }
  root.add_child("noes", noes_node);
  return root;
}
// get next power of two, from:
// https://www.geeksforgeeks.org/smallest-power-of-2-greater-than-or-equal-to-n/
const unsigned int nextPowerOf2(const unsigned int n) {
  unsigned int p = 1;
  if (n && !(n & (n - 1)))
    return n;

  while (p < n)
    p <<= 1;

  return p;
}

// rigid tumbling spectral density, assumes uniform correlation times
inline vector<MatrixXd> rigid_spectral_density(MatrixXd &dist6,
                                               vector<double> &omega,
                                               const double tau,
                                               const double total_weight) {
  vector<MatrixXd> jvec;
  // (1/4*pi) * r^(-6) * (2\tau_c/(1 + w^2 \tau_c^2 ))
  for (auto w : omega)
    jvec.push_back(dist6 * 2.0 * tau * tau /
                   ((1.0 + w * w * tau * tau) * total_weight * 4.0 * PI));
  return jvec;
}

// time conversions
const double ghz2Hz = 1e9;
const double mhz2Hz = 1e6;
const double khz2Hz = 1e3;
const double ms2s = 1e-3;
const double ns2s = 1e-9;

// second order spherical harmonic, with multiples factored out.
inline const double Y_2_0(const GCoord a, const GCoord b) {
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
  opts::WeightsOptions *wopts = new opts::WeightsOptions;
  ToolOptions *topts = new ToolOptions;

  // combine options
  opts::AggregateOptions options;
  options.add(bopts).add(sopts).add(mtopts).add(wopts).add(topts);

  // Parse the command-line.  If an error occurred, help will already
  // be displayed and it will return a FALSE value.
  if (!options.parse(argc, argv))
    exit(-1);

  // Pull the model from the options object (it will include coordinates)
  AtomicGroup model = mtopts->model;

  // Pull out the trajectory...
  pTraj traj = mtopts->trajectory;
  // if no weights file, will emit '1.0' in following code as default.
  wopts->weights->add_traj(traj);
  if (mtopts->mtraj.nframes() == 0) {
    cout << "ERROR: can't work with zero frames.\n";
    exit(-1);
  }

  // Select the desired atoms to operate over...
  AtomicGroup nuclei = selectAtoms(model, sopts->selection);
  const Eigen::Index N = (Eigen::Index)nuclei.size();
  // NOE precalculations
  const double mu0 =
      1.25663706212e4; // wikipedia, H/Angstrom (10^10 * value in H/m)
  const double gamma =
      topts->gamma * 2.0 * PI * mhz2Hz; // convert Gamma from mHz/T to Rad/s*T
  const double hbar = 1.054571817e-34;  // wikipedia, J*s
  const double N_A = 6.02214076e23;     // Wikipedia, Avogadro's Constant
  // dipolar interaction constant from Brueschweiler & Case, 1994 (eqs 33-35)
  const double xi = N_A * gamma * gamma * hbar * mu0 / (4.0 * PI);
  // this is equivalent to 0.5 * zeta from eqn 33 in B&C.
  // The 0.5 comes from factoring out of spectral density prefactors so that
  // they'd be whole numbers.
  const double zeta = xi * xi * 2.0 * PI / 10.0;
  // Larmor frequency, in Hz
  const double omega = gamma * topts->B;
  // filled with spectral densities at frqs frequencies, below.
  vector<double> frqs{0.0, omega, omega * 2.0};
  vector<MatrixXd> J;
  // matrix to hold return values of calculation in forloop
  MatrixXd sample = MatrixXd::Zero(N, N);

  // Determine which functionality the code should apply
  if (topts->corr || topts->spectral_density) {
    // get the framerate into Hertz
    const double framerate = topts->f * ghz2Hz;
    const double bin_width = topts->bin_width * khz2Hz;
    // compute fragments to Fourier Transform
    // frames_per_ft is halved to reflect causal zero padding done later.
    uint frames_per_ft, points;
    // If we need to zero pad for FFT, frames_per_ft needs to be halved for
    // desired binwidth since the following calculation will not include the
    // zeros required for causal padding.
    if (topts->corr) {
      // number of frames that would be correct for causal zero padding (2fold)
      frames_per_ft = (uint)framerate / (2 * bin_width);
      // now figure out how many points are needed for an efficient fft.
      points = nextPowerOf2(frames_per_ft);
    } else {
      frames_per_ft = (uint)framerate / bin_width;
      points = frames_per_ft; // Oscillator algorithm doesn't need zero pad
                              // since bin is centered and since rest ignored.
    }
    // Loop over multitraj to obtain individual sub-trajectory lengths;
    // the DFTs from these should be averaged separately
    vector<uint> trajlengths;
    // Create a vector of vector of frame indices for each FT.
    vector<vector<uint>> resample_FT_indices;
    auto previt = mtopts->frameList().begin();
    for (auto i = 0; i < mtopts->mtraj.size(); i++) {
      trajlengths.push_back(mtopts->mtraj.nframes(i));
      if (frames_per_ft < trajlengths.back()) {
        // taking advantage of truncation toward zero behavior of integer
        // division
        uint n_subsamples = trajlengths.back() / (frames_per_ft);
        uint bartlett_length = n_subsamples * frames_per_ft;
        uint remainder = trajlengths.back() - bartlett_length;
        for (auto j = 0; j < n_subsamples; j++) {
          // enforce a one frame overlap for each of the first 'remainder'
          // windows.
          if (j != 0 && remainder > 0) {
            remainder--;
            previt--;
          }
          vector<uint> subframes(previt, previt + j * frames_per_ft);
          resample_FT_indices.emplace_back(subframes);
          previt += j * frames_per_ft;
        }
      } else {
        cerr << "warning: trajectory " << i
             << " has fewer frames than are needed for desired bin-width.\n";
        vector<uint> trajframes(previt, previt + trajlengths.back());
        resample_FT_indices.emplace_back(trajframes);
        previt += trajlengths.back();
      }
    }

    if (topts->spectral_density) {
      vector<MatrixXd> mean_J;
      // initialize the mean J list
      for (auto i = 0; i < frqs.size(); i++) {
        MatrixXd initialized_frq = MatrixXd::Zero(N, N);
        mean_J.emplace_back(move(initialized_frq));
      }
      for (const auto subFTInds : resample_FT_indices) {
        // Magic circle oscillator precomputation:
        // Magic Circle oscillator for trackin samples like the above
        DFTMagicCircle dft(sample, frqs, framerate, subFTInds.size());
        // Now iterate over all frames in the skipped & strided trajectory
        for (const auto f : subFTInds) {
          traj->readFrame(f);
          traj->updateGroupCoords(nuclei);
          // compute second order index 0 spherical harmonics of all ij pairs.
          for (auto i = 0; i < N; i++)
            for (auto j = 0; j < i; j++)
              sample(i, j) = Y_2_0(nuclei[i]->coords(), nuclei[j]->coords());
          // feed this sample into the DFT object
          dft(sample);
        }
        // now 'finish' th
        uint zeros_count = subFTInds.size();
        while (zeros_count > 0) {
          dft(MatrixXd::Zero(N, N));
        }
        // compute spectral densities from DFT here.
        J = dft.spectral_density();
        for (auto i = 0; i < frqs.size(); i++) {
          mean_J[i] += J[i];
        }
      }
      for (auto i = 0; i < frqs.size(); i++)
        mean_J[i] /= (double) resample_FT_indices.size();
    } else {
      // create a tensor to store the FT results in
      Tensor<double, 3> mean_periodogram(points, N, N);
      // create an array to specify which dimension to FT
      Eigen::array<int, 1> dim = {(0)};
      for (const auto subFTInds : resample_FT_indices) {
        Tensor<double, 3> signal(points, N, N);
        signal.setZero();
        // Loop over the part of the multitraj corresponding to this
        // sub-transform
        for (auto f = 0; f < subFTInds.size(); f++) {
          traj->readFrame(subFTInds[f]);
          traj->updateGroupCoords(nuclei);
          // Load the tensor (3D array) with sph. harm. for this frame.
          for (auto i = 0; i < N; i++) {
            for (auto j = 0; j < i; j++) {
              signal(f, j, i) = Y_2_0(nuclei[j]->coords(), nuclei[i]->coords());
            }
          }
        }
        // accumulate periodograms (spectral densities) into the average
        auto fft_result = signal.template fft<RealPart, FFT_FORWARD>(dim);
        mean_periodogram += fft_result * fft_result;
      }
      double correction = 0;
      for (auto f = 0; f < points; f++) {
        mean_periodogram.chip(f, 0) *= 1.0 / (points - f);
      }
    }
  } else {           // use r^6 approx/averaging
    double d2 = 0.0; // stores output of squared distance
    for (auto f : mtopts->frameList()) {
      traj->readFrame(f);
      traj->updateGroupCoords(nuclei);
      // loop over ij distances
      for (auto i = 0; i < N; i++) {
        for (auto j = 0; j < i; j++) {
          d2 = nuclei[i]->coords().distance2(nuclei[j]->coords(),
                                             model.periodicBox());
          sample(i, j) += wopts->weights->get() / (d2 * d2 * d2);
        }
      }
      wopts->weights->accumulate();
    }
    // use approximate formula for spectral densitiies here.
    J = rigid_spectral_density(sample, frqs, topts->tau * ns2s,
                               wopts->weights->totalWeight());
  }
  // following van gunsteren, compute spectral density based rho and sigma
  MatrixXd rho((J[0] + (3.0 * J[1]) + (6.0 * J[2])).selfadjointView<Lower>());
  MatrixXd R(((6.0 * J[2]) - J[0]).selfadjointView<Lower>());
  R.diagonal(0) = rho.colwise().sum();
  cout << "this is R, but not in correct units:\n" << R << "\n";
  // multiply by units here
  R *= zeta;
  cout << "this is R:\n" << R << "\n";
  vector<MatrixXd> intensities;
  pt::ptree jsontree;
  jsontree.put("invocation", header);
  if (topts->isa) {
    vector<double> mix{topts->m * ms2s};
    intensities.push_back(move(R * mix[0]));
    jsontree = make_NOE_json(intensities, mix, nuclei);
  } else {
    // SelfAdjointEigensovler only looks at the lower triangle, so it is not an
    // issue that the upper triangle of the J-arithmetic above is zeros.
    SelfAdjointEigenSolver<MatrixXd> es(R);
    MatrixXd eVecs = es.eigenvectors();
    cout << "eVecs:\n" << eVecs << endl;
    MatrixXd invEvecs = eVecs.inverse();
    cout << "invEvecs:\n" << invEvecs << endl;
    cout << "R - eVecs*evals*invEvecs:\n"
         << R - (eVecs * es.eigenvalues().asDiagonal() * invEvecs) << endl;
    if (topts->M != 1.0) {
      invEvecs *= topts->M;
    }
    for (const auto time : topts->buildups) {
      intensities.push_back(eVecs *
                            (es.eigenvalues() * (-time * ms2s))
                                .array()
                                .exp()
                                .matrix()
                                .asDiagonal() *
                            invEvecs);
    }
    jsontree = make_NOE_json(intensities, topts->buildups, nuclei);
    ComputationInfo es_info = es.info();
    string comp_info_tag = "eigensolver";
    if (es_info == Success)
      jsontree.put(comp_info_tag, "Eigendecomposition successful.");
    if (es_info == NumericalIssue)
      jsontree.put(comp_info_tag,
                   "Eigendecomposition ran into a numerical issue.");
    if (es_info == NoConvergence)
      jsontree.put(comp_info_tag, "Eigendecomposition did not converge.");
    if (es_info == InvalidInput)
      jsontree.put(comp_info_tag,
                   "Eigendecomposition was given invalid input.");
  }
  pt::write_json(cout, jsontree);
}