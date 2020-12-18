
#include <iostream>
#include <regex>
#include "seekerror.h"


std::regex regStripPathFromFILELINE("\\(.*/(.*):(\\d+)\\):");

void print_exception_stack(const std::exception& e, int level)
{
    try {
        std::rethrow_if_nested(e);
    } catch(const std::exception& e2) {
        print_exception_stack(e2, level+1);
    } catch(const std::string &errStr) {
        std::cerr << "Stack_level " << level+1 << ": " << errStr << '\n';
    } catch(const char *errStr) {
        std::cerr << "Stack_level " << level+1 << ": " << errStr << '\n';
    } catch(...) {
        std::cerr << "Stack_level " << level+1 << ": Unknown exception type\n";
    }
    std::string what = regex_replace(e.what(), regStripPathFromFILELINE, "($1:$2):");
    std::cerr << "Stack_level " << level << ": " << what << '\n';
}