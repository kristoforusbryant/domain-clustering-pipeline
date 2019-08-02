#ifndef MYLIB_H_INCLUDED
#define MYLIB_H_INCLUDED


#include "containerutils.hpp"
#include "datafile.hpp"
#include "lammpsinput.hpp"
#include "lammpstraj.hpp"
#include "matrix.hpp"
#include "myclock.hpp"
#include "normalmodes.hpp"
#include "pdbio.hpp"
#include "sparsematrix.hpp"
#include "stringutils.hpp"
#include "verletlist.hpp"
#include "contactmap.hpp"

inline std::string exec(const char* cmd) {
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
    	if(fgets(buffer, 128, pipe) != NULL)
    		result += buffer;
    }
    pclose(pipe);
    return result;
}

#endif // MYLIB_H_INCLUDED
