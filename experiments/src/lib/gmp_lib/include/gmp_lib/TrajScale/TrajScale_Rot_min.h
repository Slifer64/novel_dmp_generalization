#ifndef GMP_TRAJ_ROT_MIN_SCALE_H
#define GMP_TRAJ_ROT_MIN_SCALE_H

#include <gmp_lib/TrajScale/TrajScale.h>
#include <gmp_lib/math/math.h>

namespace as64_
{

namespace gmp_
{

/** For generalizing a 3-DoF trajectory to new target/initial positions using 
 * a rotation and scaling base on
 * scale type 1 from 
 * 'A novel DMP formulation for global and frame independent spatial scaling in the task space'
 * DOI: 10.1109/RO-MAN47096.2020.9223500
 * 
 * scaling matrix: 
 * Ks = ( ||g - y0|| / ||gd - yd0|| ) * {rotation matrix that aligns gd - yd0 with g - y0 using the minimum angle of rotation}
 * 
 * where g, y0 are the new target and initial positions
 *      gd, yd0 are the demonstrated target and initial positions
 */

class TrajScale_Rot_min: public TrajScale
{
public:

  typedef std::shared_ptr<TrajScale_Rot_min> Ptr;

  /** Constructor.
   *  @param[in] n_dof: degrees of freedom.
   */
  TrajScale_Rot_min(): TrajScale(3)
  {
    this->setNominalStartFinalPos(this->Y0d, this->Ygd);
    this->setNewStartFinalPos(this->Y0, this->Yg);
  }

  enum ScaleType getScaleType() const { return TrajScale::ScaleType::ROT_MIN_SCALE; }

protected:

  arma::mat calcScaling() const
  {
    double d = arma::norm(this->Ygd - this->Y0d);
    if (d < 1e-2)
      std::cerr << "\033[1m" << "\033[33m" << "[WARNING]: " << "Scaling may be very big (on the order of "
                << std::to_string(1/d) << ") due to init and final positions being close." << "\033[0m\n";

    arma::vec nd = this->Ygd - this->Y0d;  nd = nd/arma::norm(nd);
    arma::vec n = this->Yg - this->Y0;  n = n/arma::norm(n);
    double dot_n_nd = arma::dot(n,nd);
    arma::mat R;
    if (std::fabs(std::fabs(dot_n_nd) - 1) < 1e-14)
      R = arma::mat().eye(3,3);
    else
    {
      arma::vec k = arma::cross(nd,n);
      double theta = std::acos(dot_n_nd);
      R = gmp_::axang2rotm( arma::join_vert(k, arma::vec({theta})) );
    }

    return R*( arma::norm(this->Yg - this->Y0)/arma::norm(this->Ygd - this->Y0d) );
  }

  arma::mat calcInvScaling() const
  {
    arma::vec nd = this->Ygd - this->Y0d;  nd = nd/arma::norm(nd);
    arma::vec n = this->Yg - this->Y0;  n = n/arma::norm(n);
    double dot_n_nd = arma::dot(n,nd);
    arma::mat R;
    if (std::fabs(std::fabs(dot_n_nd) - 1) < 1e-14)
      R = arma::mat().eye(3,3);
    else
    {
      arma::vec k = arma::cross(nd,n);
      double theta = std::acos(dot_n_nd);
      R = gmp_::axang2rotm( arma::join_vert(k, arma::vec({theta})) );
    }
    return R.t()*( arma::norm(this->Ygd - this->Y0d)/arma::norm(this->Yg - this->Y0) );
  }

}; // TrajScale_Rot_min

} // namespace gmp_

} // namespace as64_

#endif // GMP_TRAJ_ROT_MIN_SCALE_H
