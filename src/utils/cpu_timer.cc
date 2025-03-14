#include "cpu_timer.h"
#include <iostream>

namespace structured_fv {

ScopedCPUTimer::ScopedCPUTimer(const std::string& desc) :
  ScopedCPUTimer(desc, std::cout)
{}


ScopedCPUTimer::ScopedCPUTimer(const std::string& desc, std::ostream& os) :
  m_desc(desc),
  m_os(os)
{
  m_timer.start();
}

ScopedCPUTimer::~ScopedCPUTimer()
{
  m_timer.end();
  m_os << m_desc << ": " << m_timer.getElaspedTime() << " seconds" << std::endl;
}

}