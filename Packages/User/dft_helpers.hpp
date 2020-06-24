#ifndef DFT_HELPERS
#define DFT_HELPERS
#include <MultiTraj.hpp>
#include <cmath>
#include <exceptions.hpp>
#include <iostream> // for cerr
#include <string>
#include <vector>
namespace loos {
std::vector<std::vector<unsigned int>>
bartlett_contig(const unsigned int frames_per_ft,
                const std::vector<unsigned int> &frameList) {
  std::vector<std::vector<unsigned int>> resample_FT_indices;
  auto prevend = frameList.begin();
  auto newend = frameList.begin();
  // divvy up the whole darn list, irrespective of everything
  unsigned int n_subsamples = frameList.size() / frames_per_ft;
  if (n_subsamples == 0) {
    std::ostringstream oss;
    oss << "The pooled trajectory is not long enough for the desired "
           "bin-width\n";
    throw(LOOSError(oss.str()));
  }
  for (unsigned int i = 0; i < n_subsamples; i++) {
    newend = prevend + i * frames_per_ft;
    // std::vector<unsigned int> subFTinds(prevend, newend);
    resample_FT_indices.emplace_back(
        std::vector<unsigned int>(prevend, newend));
    prevend = newend;
  }
  return resample_FT_indices;
}

std::vector<std::vector<unsigned int>>
bartlett_resamples(const unsigned int frames_per_ft,
                   const std::vector<unsigned int> &frameList,
                   const loos::MultiTrajectory &mtraj) {
  // Loop over multitraj to obtain individual sub-trajectory lengths;
  // the DFTs from these should be averaged separately
  // Create a vector of vector of frame indices for each FT.
  std::vector<std::vector<unsigned int>> resample_FT_indices;
  // previous end of sub ft will be beginning of next ft
  auto prevend = frameList.begin();
  auto newend = frameList.begin();
  auto trajlength = 0;
  // figure out how to break up the multitraj to obtain frames_per_ft chunks
  for (uint i = 0; i < mtraj.size(); i++) {
    trajlength = mtraj.nframes(i);
    if (frames_per_ft < trajlength) {
      // Integer division truncates toward zero.
      uint n_subsamples = trajlength / frames_per_ft;
      uint total = n_subsamples * frames_per_ft;
      // compute the remainder after stretching the subsamples
      uint remainder = trajlength - total;
      // distribute remainder amongst the subsamples
      for (uint j = 0; j < n_subsamples; j++) {
        newend = prevend + j * frames_per_ft;
        // std::vector<uint> subframes(prevend, newend);
        resample_FT_indices.emplace_back(std::vector<uint>(prevend, newend));
        prevend = newend;
      }
    } else {
      std::ostringstream oss;
      oss << "warning: trajectory " << i
          << " has fewer frames than are needed for desired bin-width.\n";
      throw(LOOSError(oss.str()));
    }
  }
  return resample_FT_indices;
}
// Welch's method, meaning that windows overlap,
// multitraj is one contiguous timeseries.
std::vector<std::vector<unsigned int>>
welch_contig(const unsigned int frames_per_ft,
             const std::vector<unsigned int> &frameList,
             const float overlap_threshold = 0.2) {
  std::vector<std::vector<unsigned int>> resample_FT_indices;
  auto prevend = frameList.begin();
  auto newend = frameList.begin();
  unsigned int n_samples = static_cast<uint>(
      std::ceil(static_cast<float>(frameList.size()) / frames_per_ft));
  unsigned int overlap =
      (n_samples * frames_per_ft - frameList.size()) / n_samples;
  // complain, should upgrade to loos error.
  if (static_cast<uint>(overlap_threshold * frames_per_ft) > overlap) {
    std::ostringstream oss;
    oss.precision(3);
    oss << "overlap " << overlap << " is greater than overlap threshold ("
        << overlap_threshold << ") " << overlap_threshold << "\n";
    throw(LOOSError(oss.str()));
  }
  unsigned int remainder =
      frameList.size() - n_samples * (frames_per_ft - overlap);
  for (auto i = 0; i < n_samples; i++) {
    if (i != 0 && remainder > 0) {
      remainder--;
      newend = prevend + i * frames_per_ft - 1 - overlap;
    } else {
      newend = prevend + i * frames_per_ft - overlap;
    }
    resample_FT_indices.emplace_back(std::vector<uint>(prevend, newend));
    prevend = newend;
  }
  return resample_FT_indices;
}
// welch's method, multitraj represents independent trajectories
std::vector<std::vector<unsigned int>>
welch_resamples(const unsigned int frames_per_ft,
                const std::vector<unsigned int> &frameList,
                const MultiTrajectory &mtraj,
                const float overlap_threshold = 0.2) {
  // Loop over multitraj to obtain individual sub-trajectory lengths;
  // the DFTs from these should be averaged separately
  unsigned int trajlength;
  // Create a vector of vector of frame indices for each FT.
  std::vector<std::vector<uint>> resample_FT_indices;
  // previous end of sub ft will be beginning of next ft
  auto prevend = frameList.begin();
  auto newend = frameList.begin();
  for (uint i = 0; i < mtraj.size(); i++) {
    trajlength = mtraj.nframes(i);
    // figure out how to break up the multitraj to satisfy bin-width reqs for
    // FT.
    if (frames_per_ft < trajlength) {
      // round number of trajes that will fit at this binwidth up, then overlap
      // them
      uint n_subsamples = static_cast<uint>(
          std::ceil(static_cast<float>(trajlength) / frames_per_ft));
      uint no_overlap_total = n_subsamples * frames_per_ft;
      uint overlap = (no_overlap_total - trajlength) / n_subsamples;
      // compute the remainder after stretching the subsamples
      uint remainder = trajlength - overlap;
      // distribute remainder amongst the subsamples
      for (uint j = 0; j < n_subsamples; j++) {
        // Extra one frame overlap for each of the first 'remainder' subsamples.
        if (j != 0 && remainder > 0) {
          remainder--;
          newend = prevend + j * frames_per_ft - 1 - overlap;
        } else { // otherwise just move up by total of nonoverlappedframes.
          newend = prevend + j * frames_per_ft - overlap;
        }
        resample_FT_indices.emplace_back(
            std::vector<unsigned int>(prevend, newend));
        prevend = newend;
      }
    } else {
      throw(LOOSError("There are not enough frames in the trajectory for the "
                      "desired bin-width.\nConsider increasing bin-width."));
    }
  }
  return resample_FT_indices;
}
} // end namespace loos
#endif