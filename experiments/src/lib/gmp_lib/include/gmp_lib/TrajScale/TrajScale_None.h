#ifndef GMP_TRAJ_NONE_SCALE_H
#define GMP_TRAJ_NONE_SCALE_H

#include <gmp_lib/TrajScale/TrajScale.h>

#include <iostream>
#include <cstring>

namespace as64_
{

namespace gmp_
{

/** Produces a unit scaling (so essentially no scaling)
 */

class TrajScale_None : public TrajScale
{

// ===================================
// =======  Public Functions  ========
// ===================================
public:

  typedef std::shared_ptr<TrajScale_None> Ptr;

  /** Constructor.
   *  @param[in] n_dof: degrees of freedom.
   */
  TrajScale_None(unsigned n_dof): TrajScale(n_dof) {}

  enum ScaleType getScaleType() const { return TrajScale::ScaleType::NONE; }

protected:

  arma::mat calcScaling() const
  {
    return arma::mat().eye(n_dof, n_dof);
  }

  arma::mat calcInvScaling() const
  {
    return arma::mat().eye(n_dof, n_dof);
  }

}; // TrajScale_None

} // namespace gmp_

} // namespace as64_

#endif // GMP_TRAJ_NONE_SCALE_H
