#ifndef STRUCTURED_FINITE_VOLUME_UTILS_CPU_TIMER_H
#define STRUCTURED_FINITE_VOLUME_UTILS_CPU_TIMER_H

#include <ctime>
#include "project_defs.h"
#include <iosfwd>

namespace structured_fv {

class CPUTimer
{
  public:
  
    void start() { m_start = std::clock(); }

    double end()
    {
      m_end = std::clock();
      return getElaspedTime();
    }

    // returns time in seconds
    Real getElaspedTime() const
    {
      return static_cast<double>((m_end - m_start)) / CLOCKS_PER_SEC;
    }

  private:
    std::clock_t m_start;
    std::clock_t m_end;
};

class ScopedCPUTimer
{
  public:
    ScopedCPUTimer(const std::string& desc);

    ScopedCPUTimer(const std::string& desc, std::ostream& os);

    ~ScopedCPUTimer();

  private:
    CPUTimer m_timer;
    const std::string& m_desc;
    std::ostream& m_os;
};


}

#endif