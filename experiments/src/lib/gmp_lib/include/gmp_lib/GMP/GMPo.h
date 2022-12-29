#ifndef GMP_LIB_GMP_ORIENT_H
#define GMP_LIB_GMP_ORIENT_H

// GMPo class
//  GMP for encoding Cartesian Orientation.
//  Encodes the orientation by formulating a GMP in the quaternion log space.
//  The class provides the necessary tranformation betweeen the quaternion
//  log (and its 1st and 2nd time derivatives) and the quaternion,
//  rotational velocity and acceleration.
//

#include <gmp_lib/GMP/GMP.h>

namespace as64_
{

namespace gmp_
{

// class GMPo_Update; // forward declaration
class GMPo_IO; // forward declaration

class GMPo : public GMP
{

public:

  typedef std::shared_ptr<GMPo> Ptr;

  /** Constructs a GMP defined in the quat-log space.
   * @param[in] N_kernels: 3x1 vector with the number of kernels for each dim of the quat log.
   * @param[in] kernels_std_scaling: Scaling for std of kernels (optional, default=2).
   * \note: Each of the arguments 'N_kernels', 'D', 'K' can be scalar or a 3x1 vector.
   */
  GMPo(unsigned N_kernels=2, double kernels_std_scaling=1.0);


  /** Trains the GMPo::
   * @param[in] train_method: the training method to use, as a string ('LWR', 'LS').
   * @param[in] Time: Row vector with the timestamps of the training data points.
   * @param[in] Quat_data: Matrix with the desired unit quaternion (4x1 vector) in each column.
   * @param[out] train_error: The training error expressed as the mse error.
   */
  void train(const std::string &train_method, const arma::rowvec &Time, const arma::mat &Quat_data, arma::vec *train_error=NULL, arma::mat *Sw=NULL);


  /** Sets the initial orientation.
   * @param[in] Q0: Initial orientation (as unit quaternion).
   */
  void setQ0(const arma::vec &Q0);

  arma::vec getQ0() const { return this->Q0; }


  /** Sets the goal/target orientation.
   * @param[in] Qg: Goal/target orientation (as unit quaternion).
   */
  void setQg(const arma::vec &Qg);

  arma::vec getQd(double x) const;

  arma::vec getVd(double x, double x_dot) const;

  arma::vec getVdDot(double x, double x_dot, double x_ddot) const;

  void getRefTraj(double x, double x_dot, double x_ddot, arma::vec &Qd, arma::vec &vd, arma::vec &vd_dot) const;

  /** Returns the rotational velocity.
   * Call @update first!
   * @param[in] Q: the current orientation.
   * @return: the rotational velocity.
   */
  arma::vec getRotVel(const arma::vec &Q) const;


  /** Returns the rotational acceleration.
   * Call @update first!
   * @param[in] Q: the current orientation.
   * @param[in] yc_dot: time derivative of 'yc' coupling term (optional, default=0).
   * @return: the rotational acceleration.
   */
  arma::vec getRotAccel(const arma::vec &Q, arma::vec yc_dot={0}) const;


  /** Calculates the rotational acceleration based on the current input variables.
   * @param[in] s: Object of type @GMP_phase.
   * @param[in] Q: the current orientation.
   * @param[in] rotVel: the rotational velocity.
   * @param[in] yc: Coupling term fo 'y' state diff-equation (optional, default=0).
   * @param[in] zc: Coupling term fo 'z' state diff-equation (optional, default=0).
   * @param[in] yc_dot: time derivative of 'yc' coupling term (optional, default=0).
   * @return: the rotational acceleration.
   */
  arma::vec calcRotAccel(const gmp_::Phase &s, const arma::vec &Q, const arma::vec &rotVel, arma::vec y_c={0}, arma::vec z_c={0}, arma::vec yc_dot={0});


  /** Returns the 'y' state of the DMP based on the current orientation.
   * @param[in] Q: Current orientation (as unit quaternion).
   */
  arma::vec getY(const arma::vec &Q) const;


  /** Returns the 'z' state of the DMP based on the current rotational velocity and orientation.
   * @param[in] rotVel: Current rotational velocity.
   * @param[in] Q: Current orientation (as unit quaternion).
   */
  arma::vec getZ(const arma::vec &rotVel, const arma::vec &Q) const;

  void deepCopy(gmp_::GMPo *cp_obj) const
  {
    GMP::deepCopy(cp_obj);
    cp_obj->Q0 = this->Q0;
    cp_obj->Qd0 = this->Qd0;
  }

  // ===========================================================
  // ************      Static Public Functions      ************
  // ===========================================================

  /* Expresses a given quaternion w.r.t. the initial orientation.
   * @param[in] Q: Orientation as unit quaternion.
   * @param[in] Q0: Initial quaternion.
   * @return: The orientation w.r.t. Q0, i.e. Q1 = Q*Q0^{-1}.
   */
  static arma::vec getQ1(const arma::vec &Q, const arma::vec &Q0);


  /* Returns the log of a given orientation w.r.t. the initial orientation.
   * @param[in] Q: Orientation as unit quaternion.
   * @param[in] Q0: Initial quaternion.
   * @return: The logarithm of the Q w.r.t. Q0, i.e. q = log(Q*Q0^{-1}).
   */
  static arma::vec quat2q(const arma::vec &Q, const arma::vec &Q0);


  /* Returns the quaternion Q given the initial orientation Q0 and the log of Q w.r.t. Q0.
   * @param[in] q: Logarithm of orientation w.r.t. the initial orientation.
   * @param[in] Q0: Initial orientation.
   * @return: The orientation corresponding to log, i.e. Q = exp(q)*Q0
   */
  static arma::vec q2quat(const arma::vec &q, const arma::vec &Q0);

protected:

  // friend class GMPo_Update;
  // friend class GMPo_IO;

  friend void write(const gmp_::GMPo *gmp, gmp_::FileIO &fid, const std::string &prefix);
  friend void read(gmp_::GMPo *gmp, gmp_::FileIO &fid, const std::string &prefix);

  arma::vec Q0; // initial orientation (as unit quaternion)

  arma::vec Qd0;

private:

  static const long double zero_tol;


}; // GMPo

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_ORIENT_H
