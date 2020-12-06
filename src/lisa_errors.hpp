

#ifndef lisa_errors_hpp
#define lisa_errors_hpp


#include <fstream>
#include <stdexcept>
#include <string>
#include <string.h>


// Use '\\' for Windows
#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__) 

#define STRINGIFY(x) #x
#define TO_STRING(x) STRINGIFY(x)


///////////////////////////////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////
inline bool file_accessible(const std::string& name) {
  std::ifstream f(name);
  return f.good();
}


namespace lisa {

///////////////////////////////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////
void throw_error(const char* filename, int line, const char* message);

/// And a convenience macro.
#define THROW_ERROR(msg) throw_error(TO_STRING(__FILENAME__), __LINE__, msg)


}


#endif