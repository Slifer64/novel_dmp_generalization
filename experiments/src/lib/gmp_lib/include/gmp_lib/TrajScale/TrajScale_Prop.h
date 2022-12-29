#ifndef GMP_TRAJ_PROP_SCALE_H
#define GMP_TRAJ_PROP_SCALE_H

#include <gmp_lib/TrajScale/TrajScale.h>

#include <iostream>
#include <cstring>

namespace as64_
{

namespace gmp_
{

/** For generalizing an n-DoF trajectory to new target/initial positions using 
 * proportional scaling (as in the original DMP):
 * 
 * scaling matrix: Ks = diag( (g - y0) ./ (gd - yd0) )
 * 
 * where g, y0 are the new target and initial positions
 *      gd, yd0 are the demonstrated target and initial positions
 */

class TrajScale_Prop : public TrajScale
{

// ===================================
// =======  Public Functions  ========
// ===================================
public:

  typedef std::shared_ptr<TrajScale_Prop> Ptr;

  /** Constructor.
   *  @param[in] n_dof: degrees of freedom.
   */
  TrajScale_Prop(unsigned n_dof): TrajScale(n_dof)
  {
    this->setNominalStartFinalPos(this->Y0d, this->Ygd);
    this->setNewStartFinalPos(this->Y0, this->Yg);
  }

  enum ScaleType getScaleType() const { return TrajScale::ScaleType::PROP_SCALE; }

protected:

  arma::mat calcScaling() const
  {
    double d = arma::min( arma::abs( this->Ygd - this->Y0d ) );
    if (d < 1e-2)
      std::cerr << "\033[1m" << "\033[33m" << "[WARNING]: " << "Scaling may be very big (on the order of "
                << std::to_string(1/d) << ") due to init and final positions being close." << "\033[0m\n";

    return arma::diagmat( (this->Yg - this->Y0) / (this->Ygd - this->Y0d) );
  }

  arma::mat calcInvScaling() const
  {
    return arma::diagmat( (this->Ygd - this->Y0d) / (this->Yg - this->Y0) );
  }

}; // TrajScale_Prop

} // namespace gmp_

} // namespace as64_

#endif // GMP_TRAJ_PROP_SCALE_H
