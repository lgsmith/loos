#include <OptionsFramework.hpp>

#include <cstdlib>
#include <unistd.h>
#include <pwd.h>


namespace loos {
  namespace OptionsFramework {

    void BasicOptions::addGeneric(po::options_description& opts) {
      if (!full_help.empty())
        opts.add_options()("fullhelp", "More detailed help");

      opts.add_options()
        ("help", "Produce this message")
        ("verbosity,v", po::value<int>(&verbosity)->default_value(verbosity), "Verbosity");
    }

    std::string BasicOptions::print() const {
      std::ostringstream oss;
      oss << "verbosity=" << verbosity ;
      return(oss.str());
    }


    // Presumably, BasicOptions will be the first OptionsPackage in
    // the list of options.  This means we can catch the --fullhelp at
    // the check() stage.  If we try later (i.e. with
    // postConditions()), then another OptionsPackage may fail the
    // check (such as required args or a model name)
    bool BasicOptions::check(po::variables_map& map) {
      if (!full_help.empty()) {
        if (map.count("fullhelp")) {
          std::cout << full_help << std::endl;
          return(true);
        }
      }
      return(false);
    }


    void BasicOptions::setFullHelp(const std::string& s) {
      full_help = s;
    }

    // -------------------------------------------------------

    void OutputPrefixOptions::addGeneric(po::options_description& opts) {
      opts.add_options()
        ("prefix,p", po::value<std::string>(&prefix)->default_value(prefix), "Output prefix");
    }

    std::string OutputPrefixOptions::print() const {
      std::ostringstream oss;
      oss << "prefix='" << prefix << "'";
      return(oss.str());
    }
      
    // -------------------------------------------------------

    void BasicSelectionOptions::addGeneric(po::options_description& opts) {
      opts.add_options()
        ("selection,s", po::value<std::string>(&selection)->default_value(selection), "Which atoms to use");
    }

    std::string BasicSelectionOptions::print() const {
      std::ostringstream oss;
      oss << "selection='" << selection << "'";
      return(oss.str());
    }

    // -------------------------------------------------------

    void ModelWithCoordsOptions::addGeneric(po::options_description& opts) {
      opts.add_options()
        ("coordinates,c", po::value<std::string>(&coords_name)->default_value(coords_name), "File to use for coordinates");
    }

    void ModelWithCoordsOptions::addHidden(po::options_description& opts) {
      opts.add_options()
        ("model", po::value<std::string>(&model_name), "Model Filename");
    }


    void ModelWithCoordsOptions::addPositional(po::positional_options_description& pos) {
      pos.add("model", 1);
    }


    bool ModelWithCoordsOptions::check(po::variables_map& map) {
      return(!map.count("model"));
    }

    std::string ModelWithCoordsOptions::help() const { return("model"); }


    std::string ModelWithCoordsOptions::print() const {
      std::ostringstream oss;

      oss << boost::format("model='%s'") % model_name;
      if (!coords_name.empty())
        oss << boost::format(", coords='%s'") % coords_name;

      return(oss.str());
    }



    // -------------------------------------------------------


    void BasicTrajectoryOptions::addGeneric(po::options_description& opts) {
      opts.add_options()
        ("skip,k", po::value<unsigned int>(&skip)->default_value(skip), "Number of frames to skip")
        ("range,r", po::value<std::string>(&frame_index_spec), "Which frames to use (matlab style range)");
    };

    void BasicTrajectoryOptions::addHidden(po::options_description& opts) {
      opts.add_options()
        ("model", po::value<std::string>(&model_name), "Model filename")
        ("traj", po::value<std::string>(&traj_name), "Trajectory filename");
    }

    void BasicTrajectoryOptions::addPositional(po::positional_options_description& pos) {
      pos.add("model", 1);
      pos.add("traj", 1);
    }

    bool BasicTrajectoryOptions::check(po::variables_map& map) {
      return(! (map.count("model") && map.count("traj")));
    }

    bool BasicTrajectoryOptions::postConditions(po::variables_map& map) {
      if (skip > 0 && !frame_index_spec.empty()) {
        std::cerr << "Error- you cannot specify both a skip and a frame range...I might get confused!\n";
        return(false);
      }
      
      return(true);
    }
    
    std::string BasicTrajectoryOptions::help() const { return("model trajectory"); }
    std::string BasicTrajectoryOptions::print() const {
      std::ostringstream oss;
      oss << boost::format("model='%s', traj='%s'") % model_name % traj_name;
      if (skip > 0)
        oss << ", skip=" << skip;
      else if (!frame_index_spec.empty())
        oss << ", range='" << frame_index_spec << "'";
      
      return(oss.str());
    }
    
    // -------------------------------------------------------
    void RequiredOptions::addOption(const std::string& name, const std::string& description) {
      StringPair arg(name, description);
      if (find(arguments.begin(), arguments.end(), arg) != arguments.end()) {
        std::ostringstream oss;
        oss << "Error- duplicate command line argument requested for '" << name << "'";
        throw(LOOSError(oss.str()));
      }
      
      arguments.push_back(arg);
    }
    
    void RequiredOptions::addVariableOption(const std::string& name, const std::string& description) {
      if (vargs_set)
        throw(LOOSError("Multiple variable arguments requested"));
      
      variable_arguments = StringPair(name, description);
      vargs_set = true;
    }
    
    void RequiredOptions::addHidden(po::options_description& o) {
      for (std::vector<StringPair>::const_iterator i = arguments.begin(); i != arguments.end(); ++i)
        o.add_options()(i->first.c_str(), po::value<std::string>(), i->second.c_str());
      if (vargs_set)
        o.add_options()(variable_arguments.first.c_str(), po::value< std::vector<std::string> >(), variable_arguments.second.c_str());
    }
    
    void RequiredOptions::addPositional(po::positional_options_description& pos) {
      for (std::vector<StringPair>::const_iterator i = arguments.begin(); i != arguments.end(); ++i)
        pos.add(i->first.c_str(), 1);
      if (vargs_set)
        pos.add(variable_arguments.first.c_str(), -1);
    }
    
    
    bool RequiredOptions::check(po::variables_map& map) {
      for (std::vector<StringPair>::const_iterator i = arguments.begin(); i != arguments.end(); ++i)
        if (! map.count(i->first.c_str()) )
          return(true);
      
      if (vargs_set && !map.count(variable_arguments.first))
        return(true);
      
      return(false);
    }
    
    bool RequiredOptions::postConditions(po::variables_map& map) {
      held_map = map;
      return(true);
    }

    std::string RequiredOptions::value(const std::string& s) const {
      std::string result;
      if (held_map.count(s))
        result = held_map[s].as<std::string>();
      return(result);
    }

  
    std::vector<std::string> RequiredOptions::variableValues(const std::string& s) const {
      return(held_map[s].as< std::vector<std::string> >());
    }


    std::string RequiredOptions::help() const {
      std::string s;
      for (std::vector<StringPair>::const_iterator i = arguments.begin(); i != arguments.end(); ++i)
        s = s + " " + i->first;
      if (vargs_set)
        s = s + variable_arguments.first + " [" + variable_arguments.first + " ...]";
      return(s);
    }

    std::string RequiredOptions::print() const {
      std::ostringstream oss;
      for (std::vector<StringPair>::const_iterator i = arguments.begin(); i != arguments.end(); ++i)
        oss << i->first << "='" << held_map[i->first].as<std::string>() << "',";
      
      if (vargs_set) {
        std::vector<std::string> v = variableValues(variable_arguments.first);
        oss << variable_arguments.first << "=(";
        oss << stringVectorAsStringWithCommas(v);
        oss << ")";
      }

      return(oss.str());
    }


    // -------------------------------------------------------


    AggregateOptions& AggregateOptions::add(OptionsPackage* pack) {
      options.push_back(pack); return(*this);
    }


    void AggregateOptions::setupOptions() {
      for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
        (*i)->addGeneric(generic);

      for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
        (*i)->addHidden(hidden);

      command_line.add(generic).add(hidden);

      pos = po::positional_options_description();
      for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
        (*i)->addPositional(pos);
    }

    void AggregateOptions::showHelp() {
      std::cout << "Usage- " << program_name << " [options] ";
      for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
        std::cout << (*i)->help() << " ";
      std::cout << std::endl;
      std::cout << generic;
    }


    bool AggregateOptions::parse(int argc, char *argv[]) {
      if (program_name.empty())
        program_name = std::string(argv[0]);

      setupOptions();
      bool show_help = false;

      try {
        po::store(po::command_line_parser(argc, argv).
                  options(command_line).positional(pos).run(), vm);
        po::notify(vm);
      }
      catch (std::exception& e) {
        std::cerr << "Error- " << e.what() << std::endl;
        show_help = true;
      }

      if (!show_help)
        show_help = vm.count("help");

      if (!show_help)
        for (vOpts::iterator i = options.begin(); i != options.end() && !show_help; ++i)
          show_help = (*i)->check(vm);

      if (show_help) {
        showHelp();
        return(false);
      }

      for (vOpts::iterator i = options.begin(); i != options.end(); ++i)
        if (!(*i)->postConditions(vm)) {
          showHelp();
          return(false);
        }

      return(true);
      
    }

    std::string AggregateOptions::print() const {
      std::string result(program_name);

      result += ": ";
    
      for (vOpts::const_iterator i = options.begin(); i != options.end(); ++i) {
        result += (*i)->print();
        if (i != options.end() - 1)
          result += ", ";
      }

      return(result);
    }


    std::string AggregateOptions::header() const {
      std::ostringstream oss;

      time_t t = time(0);
      std::string timestamp(asctime(localtime(&t)));
      timestamp.erase(timestamp.length() - 1, 1);

      std::string user("UNKNOWN USER");
      struct passwd* pwd = getpwuid(getuid());
      if (pwd != 0)
        user = pwd->pw_name;

      oss << "# " << user << " (" << timestamp << ")" << std::endl;
#if defined(REVISION)
      oss << "# LOOS Version " << REVISION << std::endl;
#endif
      oss << "# " << print() << std::endl;
      return(oss.str());
    }


    std::vector<uint> assignFrameIndices(pTraj& traj, const std::string& desc, const uint skip = 0) {
      std::vector<uint> frames;

      if (desc.empty())
        for (uint i=skip; i<traj->nframes(); ++i)
          frames.push_back(i);
      else
        frames = parseRangeList<uint>(desc);

      return(frames);
          
    }

    AtomicGroup loadStructureWithCoords(const std::string model_name, const std::string coord_name = std::string("")) {
      AtomicGroup model = createSystem(model_name);
      if (!coord_name.empty()) {
        AtomicGroup coords = createSystem(coord_name);
        model.copyCoordinates(coords);
      }

      if (! model.hasCoords())
        throw(LOOSError("Error- no coordinates found in specified model(s)"));

      return(model);
    }

    std::string stringVectorAsStringWithCommas(const std::vector<std::string>& v) {
      std::ostringstream oss;
      for (uint i=0; i<v.size(); ++i)
        oss << "'" << v[i] << "'" << ((i == v.size() - 1 ) ? "" : ",");
      return(oss.str());
    }

  };
};

