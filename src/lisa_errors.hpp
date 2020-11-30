

#include <fstream>
#include <string>
#include <string.h>


// Use '\\' for Windows
#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__) 

#define STRINGIFY(x) #x
#define TO_STRING(x) STRINGIFY(x)

bool file_accessible(const std::string& name) {
  std::ifstream f(name);
  return f.good();
}



    std::string at{std::string(__FILENAME__) + ":" + std::string(TO_STRING(__LINE__))};
    throw std::runtime_error("At " + at + " file (" + source_filename + ") not accessible.\n");
