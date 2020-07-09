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
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <memory>
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
      ("threads", po::value<int>(&threads)->default_value(1),
       "Number of threads to use for FFT (used only with '-C').")
      ("tau", po::value<double>(&tau)->default_value(0),
       "set a rigid-body correlation time, in ns.")
      ("buildup-curve-range", po::value<std::string>(&buildup_range), 
       "Which mixing times (in ms) to write out for plotting (matlab style range, overrides -m).")
      ("bin-width", po::value<double>(&bin_width)->default_value(10.0),
       "Width for DFT bins, in kHz.")
      ("welch", po::bool_switch(&welch)->default_value(false),
       "Overlap subFTs such that all frames are used with selected bin-width. Otherwise, drop data from the beginning of each trajectory file.")
      ("contiguous", po::bool_switch(&contig)->default_value(false),
       "If thrown, treat each provided trajectory as part of continuous time-series. Otherwise, treat them as independent trajectories.")
      ("isa,I", po::bool_switch(&isa)->default_value(false),
       "If thrown, report relative NOEs without any spin-relaxation.")
      ("spectral-density,J", po::bool_switch(&spectral_density)->default_value(false),
       "If thrown, use spectral densities only at 0, w, and 2w for relaxation matrix elts.")
      ("correlation,C", po::bool_switch(&corr)->default_value(false),
       "if thrown, compute full DFT of result; Output spectral-density NOEs, spectra, and correlation functions.")
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
           << "You've asked for two different methods of computing \n"
              "the spectral density. Pick one.\n";
      return false;
    }
    return true;
  }
  bool isa, spectral_density, corr, welch, contig;
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
// overload for writing spectra and autocorrelations.
pt::ptree make_NOE_json(vector<MatrixXd> &intensities, vector<double> &times,
                        AtomicGroup &nuclei, Tensor<double, 3> &spectra,
                        Tensor<double, 3> &corrs, vector<double> &frq_bins,
                        vector<double> &lag_bins) {
  pt::ptree root; // we'll build this up, then return it.
  pt::ptree times_node;
  pt::ptree noes_node;
  pt::ptree spectra_node;
  pt::ptree correlations_node;
  pt::ptree frqbins_node;
  pt::ptree lagbins_node;
  // build JSON array of mixing times
  for (const auto &time : times) {
    pt::ptree time_node;
    time_node.put("", time);
    times_node.push_back(make_pair("", time_node));
  }
  root.add_child("mixing times", times_node);
  // build JSON array of frequency bins for Spectral Densities
  for (const auto &frq : frq_bins) {
    pt::ptree frq_node;
    frq_node.put("", frq);
    frqbins_node.push_back(make_pair("", frq_node));
  }
  root.add_child("frequencies", frqbins_node);
  // build JSON arry of lag-time bins for autocorrelation
  for (const auto &lag : lag_bins) {
    pt::ptree lag_node;
    lag_node.put("", lag);
    lagbins_node.push_back(move(make_pair("", lag_node)));
  }
  root.add_child("lags", lagbins_node);
  // record data for each atom, sorted by atom tag
  for (size_t i = 0; i < nuclei.size(); i++) {
    pAtom ith = nuclei[i];
    for (size_t j = 0; j < i; j++) {
      pAtom jth = nuclei[j];
      ostringstream tag;
      tag << jth->resname() << jth->resid() << jth->name() << ":"
          << ith->resname() << ith->resid() << ith->name();
      // write the spectrum and correlations to arrays tagged by the tag
      pt::ptree spectrum_node;
      pt::ptree corr_node;
      // loop over first index of tensor--this should be the bin index the
      // (i)FFT
      for (uint k = 0; k < spectra.dimension(0); k++) {
        pt::ptree spectrum_point;
        spectrum_point.put("", spectra(k, j, i));
        spectrum_node.push_back(make_pair("", spectrum_point));
        pt::ptree corr_point;
        corr_point.put("", corrs(k, j, i));
        corr_node.push_back(make_pair("", corr_point));
      }
      spectra_node.add_child(tag.str(), spectrum_node);
      correlations_node.add_child(tag.str(), corr_node);
      // writhe the correlations per lag to an array with above tag
      // loop over the first index of the tensor again
      pt::ptree noe_node;
      for (uint t = 0; t < times.size(); t++) {
        pt::ptree vol_node;
        vol_node.put("", intensities[t](j, i));
        noe_node.push_back(make_pair("", vol_node));
      }
      noes_node.add_child(tag.str(), noe_node);
    }
  }
  root.add_child("noes", noes_node);
  root.add_child("spectra", spectra_node);
  root.add_child("correlations", correlations_node);
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
  // J(w) = (1/4*pi) * <r^(-6)> * (2\tau_c/(1 + w^2 \tau_c^2 ))
  for (auto w : omega)
    jvec.emplace_back(dist6 * 2.0 * tau /
                   ((1.0 + w * w * tau * tau) * total_weight * 4.0 * PI));
  return jvec;
}

// second order spherical harmonic, with multiples factored out.
inline const double Y_2_0(const GCoord a, const GCoord b) {
  const GCoord diff = a - b;
  const double recip_r = 1.0 / diff.length();
  const double recip_r3 = recip_r * recip_r * recip_r;
  return recip_r3 * (1.5 * diff[2] * diff[2] * recip_r * recip_r - 0.5);
}

// time conversions
const double ghz2Hz = 1e9;
const double mhz2Hz = 1e6;
const double khz2Hz = 1e3;
const double ms2s = 1e-3;
const double ns2s = 1e-9;

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
  const double mu0 = 1.25663706212e-6 // wikipedia, H/m. H = T * m^2 / A
                     * 1e20           // 1e20 Angstrom/m for numerator
                     * 1e-10;         // 1e-10 m/Angstrom for denomenator
  const double gamma =
      topts->gamma * 2.0 * PI * mhz2Hz; // convert Gamma from mHz/T to Rad/(s*T)
  const double hbar = 1.054571817e-34   // wikipedia, J*s. J = kg * m^2 * Hz^2
                      * 1e20;           // * 1e20 (Angstroms/m)^2.
  // const double N_A = 6.02214076e23;     // Wikipedia, Avogadro's Constant
  // dipolar interaction constant from Brueschweiler & Case, 1994 (eqs 33-35)
  // related to equation 9, but without the factor of sqrt(12)
  const double xi = gamma * gamma * hbar * mu0 / (4.0 * PI);
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
  // make unique pointers for that may be created, depending on the following
  unique_ptr<Tensor<double, 3>> p_mean_spectra;
  unique_ptr<Tensor<double, 3>> p_corrs;
  unique_ptr<vector<double>> p_frqbins;
  unique_ptr<vector<double>> p_lagbins;
  // Determine which functionality the code should apply
  if (topts->corr || topts->spectral_density) {
    // get the framerate into Hertz
    const double fs = topts->f * ghz2Hz;
    const double bin_width = topts->bin_width * khz2Hz;
    double actual_binwidth; // store bin width after padding up to nextpwr2.
    // define these so they'll be accessible in lower scopes
    uint frames_per_ft, points;
    // If we need to zero pad for FFT, frames_per_ft needs to be halved for
    // desired binwidth since the following calculation will not include the
    // zeros required for causal padding.
    if (topts->corr) {
      // number of frames that would be correct for causal zero padding (2fold)
      frames_per_ft = (uint)fs / (2 * bin_width);
      // now figure out how many points are needed for an efficient fft.
      points = nextPowerOf2(frames_per_ft);
      actual_binwidth = fs / points;
    } else {
      frames_per_ft = (uint)fs / bin_width;
      points = frames_per_ft; // Oscillator algorithm doesn't need zero pad
                              // since bin is centered and since rest ignored.
    }
    vector<vector<uint>> resample_FT_indices;
    if (topts->welch) {
      if (topts->contig)
        resample_FT_indices = welch_contig(frames_per_ft, mtopts->frameList());
      else
        resample_FT_indices =
            welch_resamples(frames_per_ft, mtopts->frameList(), mtopts->mtraj);
    } else {
      if (topts->contig)
        resample_FT_indices =
            bartlett_contig(frames_per_ft, mtopts->frameList());
      else
        resample_FT_indices = bartlett_resamples(
            frames_per_ft, mtopts->frameList(), mtopts->mtraj);
    }
    // If we only need the three freqs from the spectral density for NOEs.
    if (topts->spectral_density) {
      // initialize the mean J list
      for (uint i = 0; i < frqs.size(); i++) {
        J.emplace_back(MatrixXd::Zero(N, N));
      }
      for (const auto subFTInds : resample_FT_indices) {
        // Magic circle oscillator precomputation:
        // Magic Circle oscillator for trackin samples like the above
        DFTMagicCircle dft(sample, frqs, fs, subFTInds.size());
        // Now iterate over all frames in the skipped & strided trajectory
        for (const auto f : subFTInds) {
          traj->readFrame(f);
          traj->updateGroupCoords(nuclei);
          // compute second order index 0 spherical harmonics of all ij pairs.
          for (auto i = 0; i < N; i++)
            for (auto j = 0; j < i; j++) {
              auto sij = Y_2_0(nuclei[i]->coords(), nuclei[j]->coords());
              sample(i, j) = sij;
            }
          // feed this sample into the DFT object
          dft(sample);
        }
        // compute spectral densities from DFT here.
        for (uint i = 0; i < frqs.size(); i++) {
          J[i] += dft.power_spectral_density()[i];
        }
      }
      for (uint i = 0; i < frqs.size(); i++)
        J[i] /= static_cast<double>(resample_FT_indices.size());
    } else { // Else we need to compute the autocorrelation, which means two
             // full FTs.
      // create a tensor to store the FT results in
      Tensor<double, 3> mean_periodogram(points, N, N);
      mean_periodogram.setZero();
      // have the outer-scope unique pointer point at it, for saving data l8r
      p_mean_spectra.reset(&mean_periodogram);
      // create an array to specify which dimension to FT
      Eigen::array<int, 1> dim = {(0)};
      for (const auto subFTInds : resample_FT_indices) {
        Tensor<double, 3> signal(points, N, N);

        signal.setZero();
        // Loop over the part of the multitraj corresponding to this
        // sub-transform
        for (uint f = 0; f < subFTInds.size(); f++) {
          traj->readFrame(subFTInds[f]);
          traj->updateGroupCoords(nuclei);
          // Load the tensor (3D array) with sph. harm. for this frame.
          for (auto i = N - 1; i != 0; i--) {
            for (auto j = 0; j < i; j++) {
              signal(f, i, j) = Y_2_0(nuclei[i]->coords(), nuclei[j]->coords());
            }
          }
        }
        // accumulate periodograms (spectral densities) into the average
        mean_periodogram +=
            (signal.template fft<RealPart, FFT_FORWARD>(dim)).abs().square();
      }
      // define a scale term equal to the timestep sq / length_of_traj
      // making the periodogram a Power Spectral Density
      double scale = 1.0 / (static_cast<double>(points) * fs * static_cast<double>(resample_FT_indices.size()));     
      // scale and normalize the mean_periodogram
      mean_periodogram = mean_periodogram * scale;
      
      // figure out which bins are needed to fill out J, above.
      for (const auto w : frqs) {
        // casting to unsigned int truncates, instead of rounding
        uint k = static_cast<unsigned int>(w / (2 * PI * fs));
        MatrixXd amplitude = MatrixXd::Zero(N, N);
        // cool kids would do this with a map to the chip...
        for (auto i = N - 1; i != 0; i--)
          for (auto j = 0; j < i+1; j++)
            amplitude(i, j) = mean_periodogram(k, i, j);
        J.emplace_back(move(amplitude));
      }
      for (const auto j : J)
        cout << "An amplitude record:\n" << j << endl;

      // finish calculating the discrete correlation
      Tensor<double, 3> correlation =
          mean_periodogram.template fft<RealPart, FFT_REVERSE>(dim);
      // build the frq and lag bin vectors for plotting
      vector<double> lag_bins;
      vector<double> frq_bins;
      for (uint k = 0; k < points; k++) {
        // 
        correlation.chip(k, 0) =
            correlation.chip(k, 0) * 1.0 / static_cast<double>(points - k);
        // lags should be distributed as reciprocal of bin-width
        lag_bins.push_back(k / actual_binwidth);
        // lower bin edges should be a function of bin-width
        frq_bins.push_back(k * actual_binwidth);
      }
      // point p_corrs at correlation, for writing data later.
      p_corrs.reset(&correlation);
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
    // use approximate formula for spectral densities here.
    J = rigid_spectral_density(sample, frqs, topts->tau * ns2s,
                               wopts->weights->totalWeight());
  }
  // following van Gunsteren, compute spectral density based rho and sigma
  MatrixXd rho((J[0] + (3.0 * J[1]) + (6.0 * J[2])).selfadjointView<Lower>());
  MatrixXd R(((6.0 * J[2]) - J[0]).selfadjointView<Lower>());
  R.diagonal(0) = rho.colwise().sum();
  cerr << "this is R, but not in correct units:\n" << R << "\n";
  // multiply by units here
  R *= zeta;
  cerr << "this is R:\n" << R << "\n";
  vector<MatrixXd> intensities;
  pt::ptree jsontree;
  if (topts->isa) {
    vector<double> mix{topts->m};
    intensities.push_back(move(R * ms2s * mix[0]));
    jsontree = make_NOE_json(intensities, mix, nuclei);
  } else {
    // SelfAdjointEigensovler only looks at the lower triangle, so it is not an
    // issue that the upper triangle of the J-arithmetic above is zeros.
    SelfAdjointEigenSolver<MatrixXd> es(R);
    MatrixXd eVecs = es.eigenvectors();
    cerr << "eVecs:\n" << eVecs << endl;
    MatrixXd invEvecs = eVecs.inverse();
    cerr << "invEvecs:\n" << invEvecs << endl;
    cerr << "eigenvalues:\n" << es.eigenvalues() << endl;
    cerr << "R - eVecs*evals*invEvecs:\n"
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
    if (topts->corr)
      jsontree =
          make_NOE_json(intensities, topts->buildups, nuclei, *p_mean_spectra,
                        *p_corrs, *p_frqbins, *p_lagbins);
    else
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
  jsontree.put("invocation", header);
  pt::write_json(cout, jsontree);
}