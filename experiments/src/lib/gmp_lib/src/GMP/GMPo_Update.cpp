#include <gmp_lib/GMP/GMPo_Update.h>

// N-DoF GMP Update class
//

#include <gmp_lib/GMP/GMPo.h>
#include <gmp_lib/GMP/GMP_Update.h>
#include <gmp_lib/math/quaternions.h>

namespace as64_
{

namespace gmp_
{

  GMPo_Update::GMPo_Update(GMPo *gmp): GMP_Update(gmp)
  {
    this->gmp_o = gmp;
  }

  void GMPo_Update::updateQuat(double x, const arma::vec &Q, double r_n)
  {
    arma::vec y = GMPo::quat2q(Q, this->gmp_o->getQ0());
    updatePos(x, y, r_n);
  }

  void GMPo_Update::updateRotVel(double x, double x_dot, const arma::vec &rotVel, const arma::vec &Q, double r_n)
  {
    arma::vec Q1 = GMPo::getQ1(Q, this->gmp_o->getQ0());
    arma::vec y_dot = gmp_::rotVel_to_qLogDot(rotVel, Q1);
    updateVel(x, x_dot, y_dot, r_n);
  }

  void GMPo_Update::updateRotAccel(double x, double x_dot, double x_ddot, const arma::vec &rotAccel, const arma::vec &rotVel, const arma::vec &Q, double r_n)
  {
    arma::vec Q1 = GMPo::getQ1(Q, this->gmp_o->getQ0());
    arma::vec y_ddot = gmp_::rotAccel_to_qLogDDot(rotAccel, rotVel, Q1);
    updateAccel(x, x_dot, x_ddot, y_ddot, r_n);
  }

} // namespace gmp_

} // namespace as64_
