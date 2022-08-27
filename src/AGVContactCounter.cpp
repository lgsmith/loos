#include <AGVContactCounter.hpp>

namespace loos {

  std::vector<uint> AGVContactCounter::operator()(void) {
    // update contactor coords
    for (size_t i=0; i<contactors.size(); i++)
      contactor_coords[i] = contactors[i]->coords();
    // count contacts
    for (size_t i=0; i<contact_count_list.size(); i++){
      contact_count_list[i] = contacted[i].totalContacts(dist, contactor_coords);
    }
    return(contact_count_list);
  }
  std::vector<uint> AGVContactCounter::operator()(GCoord& box) {
    // update contactor coords
    for (size_t i=0; i<contactors.size(); i++)
      contactor_coords[i] = contactors[i]->coords();
    for (size_t i=0; i<contact_count_list.size(); i++){
      contact_count_list[i] = contacted[i].totalContacts(dist, contactor_coords, box);
    }
    return(contact_count_list);
  }  
} // namespace loos