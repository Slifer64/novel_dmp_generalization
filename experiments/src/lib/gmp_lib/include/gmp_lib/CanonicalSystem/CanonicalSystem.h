
#ifndef AS64_GMP_CANONICAL_SYSTEM_H
#define AS64_GMP_CANONICAL_SYSTEM_H

namespace as64_
{

namespace gmp_
{


class CanonicalSystem
{
    
// ============================
// ========  public  ==========
// ============================

public:

  double s; ///< phase
  double s_dot; ///< phase 1st time derivative
  double sd_dot; ///< desired speed of the phase evolution
    
  /** Constructor
   * @parram[in] T: time duration for the phase to go from 0 to 1.
    * @parram[in] Ds: rate of tracking the time duration. (default = 30).
    *              Set to inf to allow instant (discontinuous changes)
    */
  CanonicalSystem(double T, double Ds=30);
    
  /** Integrate the phase.
   * @parram[in] t0: initial time instant.
   * @parram[in] tf: final time instant.
   */
  void integrate(double t0, double tf);
    
  /** Set the time duration for the phase to reach 1, given the current time instant.
   * @parram[in] T: time duration
   * @parram[in] t: current time instant (default = 0).
   */
  void setDuration(double T, double t=0);
  
  /** Set the remaing time duration for the phase to reach 1.
   * @parram[in] T: reaming time duration
   */
  void setRemainingDuration(double T);
      
  /** Set the duration
   * @parram[in] s: phase.
   * @parram[in] s_dot: phase 1st time derivative.
   * @return: the 2nd time derivative of the phase
   */
  double getPhaseDDot(double s, double s_dot);
  
  /** Set state values to their initial values.
   */ 
  void reset();

// ===============================
// ========  protected  ==========
// ===============================

protected:

  double s0; ///< initial value of the phase
  double sf; ///< final value of the phase

  double Ds; ///< rate of tracking the desired speed of the phase evolution
    
}; // class CanonicalSystem

} // namespace gmp_

    
} // namespace as64_

#endif // AS64_GMP_CANONICAL_SYSTEM_H