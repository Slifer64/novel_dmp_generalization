#ifndef GMP_TRAJ_ROT_WB_SCALE_H
#define GMP_TRAJ_ROT_WB_SCALE_H

#include <gmp_lib/TrajScale/TrajScale.h>
#include <gmp_lib/math/math.h>

namespace as64_
{

namespace gmp_
{

/** For generalizing a 3-DoF trajectory to new target/initial positions using 
 * a rotation and scaling base on
 * scale type 3 from 
 * 'A novel DMP formulation for global and frame independent spatial scaling in the task space'
 * DOI: 10.1109/RO-MAN47096.2020.9223500
 * 
 * scaling matrix: 
 * Ks = ( ||g - y0|| / ||gd - yd0|| ) * {rotation matrix that aligns gd - yd0 with g - y0 and the rotation axis is the closest to the workbench normal}
 * 
 * where g, y0 are the new target and initial positions
 *      gd, yd0 are the demonstrated target and initial positions
 */

class TrajScale_Rot_wb : public TrajScale
{
public:

  typedef std::shared_ptr<TrajScale_Rot_wb> Ptr;

  /** Constructor.
   *  @param[in] n_dof: degrees of freedom.
   */
  TrajScale_Rot_wb(): TrajScale(3)
  {
    this->setNominalStartFinalPos(this->Y0d, this->Ygd);
    this->setNewStartFinalPos(this->Y0, this->Yg);
    this->setWorkBenchNormal({0,0,1});
  }

  /** Sets the normal to the work bench vector. Applies only to 'ROT_WB_SCALE'.
   *  @param[in] n_wb: 3x1 vector of the normal to the workbench.
   */
  void setWorkBenchNormal(const arma::vec &n_wb)
  {
    this->n_wb = n_wb;

    this->calcScaling();
    this->calcInvScaling();
  }

  arma::vec getWorkBenchNormal() const { return n_wb; }

  ScaleType getScaleType() const { return TrajScale::ScaleType::ROT_WB_SCALE; }

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
    if (std::fabs(std::fabs(dot_n_nd) - 1) < 1e-12)
      R = arma::mat().eye(3,3);
    else
    {
      arma::vec p = (n - nd) / arma::norm(n - nd);
      arma::vec k = this->n_wb - arma::dot(p,this->n_wb)*p;
      k = k / arma::norm(k);

      arma::vec nd2 = nd - arma::dot(k,nd)*k;
      arma::vec n2 = n - arma::dot(k,n)*k;
      double cos_n = arma::dot(n2,nd2) / (arma::norm(n2)*arma::norm(nd2));
      double theta;
      if ( std::fabs(std::fabs(cos_n)-1) < 1e-15 ) theta = arma::datum::pi;
      else theta = std::acos( cos_n );

      R = gmp_::axang2rotm( arma::join_vert(k, arma::vec({theta})) );

      if (arma::norm(n - R*nd) > 1e-8) R = gmp_::axang2rotm( arma::join_vert(k, arma::vec({-theta})) );
    }

    if (arma::norm(n - R*nd) > 1e-6) throw std::runtime_error("[TrajScale_Rot_wb::calcScaling]: R produces significant error...") ;

    return R * ( arma::norm(this->Yg - this->Y0)/arma::norm(this->Ygd - this->Y0d) );

  }

  arma::mat calcInvScaling() const
  {
    return arma::inv(this->calcScaling());
  }


private:
  arma::vec n_wb; ///< workbench normal

}; // TrajScale_Rot_wb

} // namespace gmp_

} // namespace as64_

#endif // GMP_TRAJ_ROT_WB_SCALE_H
