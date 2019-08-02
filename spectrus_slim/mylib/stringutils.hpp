/********* stringutils.hpp **************
Useful string functions

Author: Guido Polles
Date: 30/05/2014
*************************************/

#ifndef __STRINGUTILS_HPP__
#define __STRINGUTILS_HPP__

#define MAX_LINE_LEN 10000

#include <iostream>
#include <vector>
#include <cstring>
#include <cstdio>
#include <string>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>

#include <sstream>

namespace mylib{

inline bool is_empty_line(const char* line){
  static const char *emptyline_detector = " \t\n";
  return strspn(line, emptyline_detector) == strlen(line);
}

inline size_t strtrim(char *out, size_t len, const char *str)
{
  if(len == 0)
    return 0;

  const char *end;
  size_t out_size;

  // Trim leading space
  while(isspace(*str)) str++;

  if(*str == 0)  // All spaces?
  {
    *out = 0;
    return 1;
  }

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;
  end++;

  // Set output size to minimum of trimmed string length and buffer size minus 1
  out_size = (end - str) < len-1 ? (end - str) : len-1;

  // Copy trimmed string and add null terminator
  memcpy(out, str, out_size);
  out[out_size] = 0;

  return out_size;
}

inline std::vector<std::string> strsplit(char* str, const char* delim = " \n"){
    std::vector<std::string> v;
    char *pch;
    pch = strtok(str,delim);
    while (pch != NULL){
        v.push_back(std::string(pch));
        pch = strtok (NULL, delim);
    }
    return v;
}


inline int atoi(const char* input){
    int temp;
    std::istringstream ssout(input);
    ssout>>temp;
    return temp;
}
inline int atoi(std::string input){
    int temp;
    std::istringstream ssout(input);
    ssout>>temp;
    return temp;
}
inline double atof(const char* input){
    double temp;
    std::istringstream ssout(input);
    ssout>>temp;
    return temp;
}
inline double atof(std::string input){
    double temp;
    std::istringstream ssout(input);
    ssout>>temp;
    return temp;
}

inline std::string itoa(long int input){
    std::stringstream ss;
    ss<<input;
    return ss.str();
}

inline std::string ftoa(double input){
    std::stringstream ss;
    ss<<input;
    return ss.str();
}



// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}

} // namespace mylib
#endif
