#ifndef __MYCLOCK_HPP__
#define __MYCLOCK_HPP__

#include <ctime>
#include <vector>
#include <string>
#include <iomanip>

#ifndef NULL
#define NULL 0
#endif

namespace mylib{


  struct clk_t{
    time_t start;
    time_t end;
    bool running;
    std::string label;
  };

  class Clock {
  public:
    Clock(){}
    virtual ~Clock(){}

   void start(const std::string& label) {
    int idx = find(label);
    if( idx != -1){
      _clocks[idx].running = true;
      _clocks[idx].start = time(NULL);
    }else{
      clk_t clk;
      clk.running = true;
      clk.label = label;
      clk.end = 0;
      _clocks.push_back(clk);
      _clocks.back().start = time(NULL);
    }
  }

   void stop(const std::string& label) {
    int idx = find(label);
    if( idx != -1){
      _clocks[idx].end = time(NULL);
      _clocks[idx].running = false;
    }else{
      std::cerr << "Clock::stop() error. No clock labeled \""
          << label << "\"" << std::endl;
    }
  }

  void print(const std::string& label, std::ostream& stream = std::cout) {
    int idx = find(label);

    if( idx != -1){
      long int elapsed = _clocks[idx].running ? time(NULL) - _clocks[idx].start:
          _clocks[idx].end - _clocks[idx].start;
      int seconds = elapsed % 60;
      elapsed/=60;
      int minutes = elapsed % 60;
      elapsed/=60;
      int hours = elapsed % 24;
      elapsed/=24;
      int days = elapsed;
      stream << "Timing for " << _clocks[idx].label << ": "
          << std::setw(3) << days << " d "
          << std::setw(3) << hours << " h "
          << std::setw(3) << minutes << " m "
          << std::setw(3) << seconds << " s " << std::endl;
    }else{
      std::cerr << "Clock::print() error. No clock labeled \""
          << label << "\"" << std::endl;
    }
  }

  long int getTime(const std::string& label) {
    int idx = find(label);
    if( idx != -1){
      return _clocks[idx].running ? time(NULL) - _clocks[idx].start:
          _clocks[idx].end - _clocks[idx].start;
    }else{
      std::cerr << "Clock::getTime() error. No clock labeled \""
          << label << "\"" << std::endl;
      return -1;
    }
  }

  private:
  int find(const std::string& label) const {
    int index = -1;
    for (int i = 0; i < int(_clocks.size()); ++i) {
      if(label==_clocks[i].label){
        index = i;
        break;
      }
    }
    return index;
  }



  private:
    std::vector<clk_t> _clocks;
  };



  class MyChrono{
  public:
    void start(){startc=clock();}
    void stop(){stopc=clock();}
    clock_t elapsed(){return stopc-startc;}

  private:
    clock_t startc;
    clock_t stopc;
  };
} // namespace mylib
#endif
