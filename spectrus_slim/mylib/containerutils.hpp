#ifndef __CONTAINERUTILS_HPP__
#define __CONTAINERUTILS_HPP__

#include <vector>
#include <list>
#include <algorithm>

namespace mylib{

template<class T>
inline bool is_in(const T& item,const std::vector<T>& vector)
{ return std::find(vector.begin(), vector.end(), item)!=vector.end(); }

template<class T>
inline bool is_in(const T& item,const std::list<T>& list)
{ return std::find(list.begin(), list.end(), item)!=list.end(); }

template<class T>
T max(const T& item, const std::vector<T>& vector)
{
  T max = vector[0];
  for (size_t i = 1; i < vector.size(); ++i){
    if (vector[i]>max) max=vector[i];
  }
  return max;
}

template<class T>
T max(const T& item,const std::list<T>& vector)
{
  typename std::list<T>::iterator it = vector.begin();
  T max = *it;
  for ( ++it; it != vector.end() ; ++it){
    if (*it>max) max=*it;
  }
  return max;
}

}// namespace mylib

#endif
