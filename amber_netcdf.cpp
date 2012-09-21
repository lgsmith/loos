// (c) 2012 Tod D. Romo, Grossfield Lab, URMC

#include <amber_netcdf.hpp>
#include <AtomicGroup.hpp>

namespace loos {

  void AmberNetcdf::init(const char* name, const int natoms) {
    int retval;

    retval = nc_open(name, NC_NOWRITE, &_ncid);
    if (retval)
      throw(AmberNetcdfOpenError());

    readGlobalAttributes();

    // Verify # of atoms match...
    int atom_id;
    retval = nc_inq_dimid(_ncid, "atom", &atom_id);
    if (retval)
      throw(AmberNetcdfError("Error reading atom id", retval));
    retval = nc_inq_dimlen(_ncid, atom_id, &_natoms);
    if (retval)
      throw(AmberNetcdfError("Error reading atom length", retval));
    if (_natoms != natoms)
      throw(LOOSError("AmberNetcdf has different number of atoms than expected"));


    // Get nframes
    int frame_id;
    retval = nc_inq_dimid(_ncid, "frame", frame_id);
    if (retval)
      throw(AmberNetcdfError("Error reading frame information", retval));
    retval = nc_inq_dimlen(_ncid, frame_id, &_nframes);

    // Check for periodic cells...
    retval = nc_inq_varid(_ncid, "cell_lengths", &_cell_lengths_id);
    if (!retval) {
      _periodic = true;
      retval = nc_inq_vartype(_ncid, _cell_lengths_id, &_box_type);
      if (retval)
        throw(AmberNetcdfError("Error reading periodic cell data type", retval));
      if (_box_data)
        throw(LOOSError("AmberNetcdf::init() internal error - box_data already exists"));
      switch(_box_type) {
      case NC_FLOAT: _box_data = new float[3]; break;
      case NC_DOUBLE: _box_data = new double[3]; break;
      default: throw(AmberNetcdfError("Only float and double supported for cell_lengths variable"));
      }
    }

    // Any additional validation should go here...

    // Get coord-id for later use...
    retval = nc_inq_varid(_ncid, "coordinates", &_coord_id);
    if (retval)
      throw(AmberNetcdfError("Error getting id for coordinates", retval));
    retval = nc_inq_vartype(_ncid, _coord_id, &_coord_type);
    if (retval)
      throw(AmberNetcdfError("Error getting data type for coordinates", retval));
    switch(_coord_type) {
    case NC_FLOAT: _coord_data = new float[_natoms * 3]; break;
    case NC_DOUBLE: _coord_data = new double[_natoms * 3]; break;
    default: throw(AmberNetcdfError("Only float and double supported for coordinates")); break;
    }
    
  }


  void AmberNetcdf::readRawFrame(const uint frameno) {
    size_t start[3] = { 0, 0, 0 };
    size_t count[3] = {1, 1, 3};


    // Read coordinates first...
    start[0] = frameno;
    count[1] = _natoms * 3;

    int retval;
    switch(_coord_type) {
    case NC_FLOAT: retval = nc_get_vara_float(_ncid, _coord_id, start, count, _coord_data); break;
    case NC_DOUBLE: retval = nc_get_vara_double(_ncid, _coord_id, start, count, _coord_data); break;
    default: throw(AmberNetcdfError("Only float and double data types supported for coordinates."));
    }

    if (retval)
      throw(AmberNetcdfError("Error while reading Amber netcdf frame", retval));

    // Now get box if present...
    if (_periodic) {
      count[1] = 3;
      switch(_box_type) {
      case NC_FLOAT: retval = nc_get_vara_float(_ncid, _cell_lengths_id, start, count, _box_data); break;
      case NC_DOUBLE: retval = nc_get_vara_double(_ncid, _cell_lengths_id, start, count, _box_data); break;
      default: throw(AmberNetcdfError("Only float and double data type supported for periodic boxes."));
      }
      ir (retval)
        throw(AmberNetcdfError("Error while reading Amber netcdf periodic box", retval));
      
    }
    
  }


  void AmberNetcdf::seekNextFrameImpl() {
    ++_current_frame;
  }

  void AmberNetcdf::seekFrameImpl(const uint i) {
    _current_frame = i;
  }

  bool AmberNetcdf::parseFrame() {
    if (_current_frame >= _nframes)
      return(false);

    readRawFrame(_current_frame);
  }

  void AmberNetcdf::updateGroupCoords(AtomicGroup& g) {
    if (g.size() != _natoms)
      throw(AmberNetcdfError("Cannot use AmberNetcdf::updateGroupCoords() when passed group has different number of atoms"));

    uint j=0;
    for (uint i=0; i<_natoms; ++i, j += 3)
      g[i]->coords(GCoord(_coord_data[j], _coord_data[j+1], _coord_data[j+2]));

    if (_periodic)
      g.periodicBox(GCoord(_box_data[0], _box_data[1], _box_data[2]));
  }


  void AmberNetcdf::readGlobalAttributes() {
    int retval;

    _title = readGlobalAttribute("title");
    _application = readGlobalAttribute("application");
    _program = readGlobalAttribute("program");
    _programVersion = readGlobalAttribute("programVersion");
    _conventions = readGlobalAttribute("Conventions");
    _conventionVersion = readGlobalAttribute("ConventionVersion");
  }


  // Will return an emptry string if the attribute is not found
  string AmberNetcdf::readGlobalAttribute(const std::string& name) {
    size_t len;
    
    retval = nc_inq_attname(_ncid, NC_GLOBAL, name.c_str(), &len);
    if (retval)
      return(string());

    nc_type type;
    retval = nc_inq_atttype(_ncid, NC_GLOBAL, name.c_str(), &type);
    if (type != NC_CHAR)
      throw(AmberNetcdfTypeError("Only character data is supported for global attributes"));
    

    char* buf = new char[len+1];
    retval = nc_get_att_text(_ncid, NC_GLOBAL, name.c_str(), buf);
    if (retval) {
      delete[] buf;
      throw(AmberNetcdfError("Error reading attribute " + name));
    }

    return(string(buf));
  }





};