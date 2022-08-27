#if !defined(LOOS_AGVCONTACT_COUNTER_HPP)
#define LOOS_AGVCONTACT_COUNTER_HPP

#include <AtomicGroup.hpp>
#include <loos_defs.hpp>
#include <vector>


namespace loos {
class AGVContactCounter {
public:
  AGVContactCounter(std::vector<AtomicGroup>& contacted, AtomicGroup& contactors, const greal dist) : contacted(contacted), contactors(contactors), dist(dist), contactor_coords(contactors.size(), GCoord()), contact_count_list(contacted.size(), 0) {}

  std::vector<uint> operator()(void);
  std::vector<uint> operator()(GCoord& box);
private:
  std::vector<AtomicGroup> contacted;
  AtomicGroup contactors;
  const greal dist; 

  std::vector<GCoord> contactor_coords;
  std::vector<uint> contact_count_list;
  
};

} // namespace loos
#endif