

#include <sstream>


namespace lisa {


void throw_error(const char* filename, int line, const char* message) {

  std::stringstream ss;
  ss << "Error at " << filename << ":" << line << " " << message;

  throw std::runtime_error(ss.str());  

}


}
