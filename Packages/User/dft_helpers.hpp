#ifndef DFT_HELPERS
#define DFT_HELPERS
#include <vector>

// a function for computing the details that need to go into
// takes a desired bin width, a frame rate, and a traj_length
// returns a vector of index vectors to be used for each sub-transform
const unsigned int binwidth2welch_offset(const double desired_binwidth,
                                         const double framerate,
                                         const double resample_overlap) {
  const unsigned int frames_per_ft = (unsigned int)framerate / desired_binwidth;
  const unsigned int overlap_offset =
      (unsigned int)frames_per_ft * resample_overlap;

  return overlap_offset;
}

// std::vector<std::vector<unsigned int>>
// build_framelists(const unsigned int offset, const std::vector<unsigned int>) {
//   std::vector<std::vector<unsigned int>> framelists;
//   for (auto i = 0; i < offset; i++){
    
//   }
// }

#endif