#ifndef GMP_LIB_GMP_ORIENT_UPDATE_H
#define GMP_LIB_GMP_ORIENT_UPDATE_H

// N-DoF GMP Update class
//

#include <gmp_lib/GMP/GMPo.h>
#include <gmp_lib/GMP/GMP_Update.h>
#include <gmp_lib/math/quaternions.h>

namespace as64_
{

namespace gmp_
{

class GMPo_Update : public GMP_Update
{
public:

  typedef std::shared_ptr<GMPo_Update> Ptr;

  // GMP constructor.
  //  @param[in] gmp: n_DoF dmp.
  GMPo_Update(GMPo *gmp);

  // function initSigmaw(this)

  // function initSigmaWfromMsr(this, x_data)

  // function enableSigmawUpdate(this, flag)

  // function setSigmaW(this, Sw)

  // function Sw = getSigmaW(this)

  // function setMsrNoiseVar(this, rv)


  // ==================================================
  // ===============   Online update  =================
  // ==================================================

  // function updatePos(this, x, y, r_n)

  // function updateVel(this, x, x_dot, y_dot, r_n)

  // function updateAccel(this, x, x_dot, x_ddot, y_ddot, r_n)

  void updateQuat(double x, const arma::vec &Q, double r_n=nan(""));

  void updateRotVel(double x, double x_dot, const arma::vec &rotVel, const arma::vec &Q, double r_n=nan(""));

  void updateRotAccel(double x, double x_dot, double x_ddot, const arma::vec &rotAccel, const arma::vec &rotVel, const arma::vec &Q, double r_n=nan(""));

  // function updateWeights(this, s, Z, type, r_n)

  // function plotWeightsCovariance(this, ax)

private:

  GMPo *gmp_o;

}; // GMP_Update

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_ORIENT_UPDATE_H
