
#include <iostream>
#include <regex>
#include "seekerror.h"

// Regex to remove path from the __FILE__ string
std::regex regStripPathFromFILELINE("\\(.*/(.*):(\\d+)\\):");

std::string print_exception_stack(const std::exception& e, int level)
{
    std::string trace;
    try {
        std::rethrow_if_nested(e);
    } catch(const std::exception& e2) {
        std::string ctrace = print_exception_stack(e2, level+1);
        trace += ctrace + " | ";
    } catch(const std::string &errStr) {
        std::string ctrace = "Stack_level " + std::to_string(level+1) + ": " + errStr;
        std::cerr << ctrace << std::endl;
        trace += ctrace + " | ";
    } catch(const char *errStr) {
        std::string ctrace = "Stack_level " + std::to_string(level+1) + ": " + errStr;
        std::cerr << ctrace << std::endl;
        trace += ctrace + " | ";
    } catch(...) {
        std::string ctrace = "Stack_level " + std::to_string(level+1) + ": Unknown exception type";
        std::cerr << ctrace << std::endl;
        trace += ctrace + " | ";
    }
    std::string what = regex_replace(e.what(), regStripPathFromFILELINE, "($1:$2):");
    std::string ttrace = "Stack_level " + std::to_string(level) + ": " + what;
    std::cerr << ttrace << std::endl;
    trace += ttrace;
    return trace;
}