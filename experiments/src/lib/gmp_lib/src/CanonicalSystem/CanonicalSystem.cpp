#include <gmp_lib/CanonicalSystem/CanonicalSystem.h>

#include <cstdlib>
#include <cmath>
#include <array>
#include <exception>

#include <boost/numeric/odeint.hpp>

namespace as64_
{

namespace gmp_
{

typedef std::array<double,2> can_sys_state;
typedef boost::numeric::odeint::runge_kutta_cash_karp54< can_sys_state > can_sys_error_stepper_type;
typedef boost::numeric::odeint::controlled_runge_kutta< can_sys_error_stepper_type > can_sys_controlled_stepper_type;

CanonicalSystem::CanonicalSystem(double T, double Ds)
{    
  this->s0 = 0;
  this->sf = 1;
  this->Ds = Ds;

  this->sd_dot = 1/T;
    
  this->s = 0;
  this->s_dot = this->sd_dot;
}
  

void CanonicalSystem::integrate(double t0, double tf)
{
  static can_sys_controlled_stepper_type ctrl_stepper;

  if (tf < t0) throw std::runtime_error("Cannot integrate backwards in time: t0 > tf");

  if (std::isinf(this->Ds))
  {
    this->s_dot = this->sd_dot;
    this->s = this->s + this->s_dot*(tf - t0);
    
    if (this->s > this->sf) this->s = this->sf; 
    else if (this->s < this->s0) this->s = this->s0;

    return;
  }

  auto ode_fun = [this](const can_sys_state &state, can_sys_state &state_dot, double t)
  { 
    state_dot[0] = state[1];
    state_dot[1] = this->getPhaseDDot(state[0], state[1]); 
  };
  
  can_sys_state state = {s, s_dot};
  boost::numeric::odeint::integrate_adaptive(ctrl_stepper, ode_fun, state, t0, tf, tf-t0);
  s = state[0];
  s_dot = state[1];
}
  
void CanonicalSystem::setDuration(double T, double t)
{
  if (T < t) throw std::runtime_error("The current time has already exceeded the duration");
      
  this->sd_dot = (this->sf - this->s) / (T - t);
}

void CanonicalSystem::setRemainingDuration(double T)
{      
  this->sd_dot = (this->sf - this->s) / T;
}

double CanonicalSystem::getPhaseDDot(double s, double s_dot)
{
  double s_ddot = 0;
  if (s>=this->s0 && s < this->sf) s_ddot = -this->Ds*(s_dot - this->sd_dot);
  else if (s>=1)                   s_ddot = -400*s_dot - 1000*(s-this->sf);
  else                             s_ddot = -400*s_dot - 1000*(s-this->s0);
  
  return s_ddot;
}

void CanonicalSystem::reset()
{
  this->s = this->s0;
  this->s_dot = this->sd_dot;
}


} // namespace gmp_

    
} // namespace as64_