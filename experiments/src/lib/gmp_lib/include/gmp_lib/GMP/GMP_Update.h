#ifndef GMP_LIB_GMP_UPDATE_H
#define GMP_LIB_GMP_UPDATE_H

// N-DoF GMP Update class
//

#include <gmp_lib/GMP/GMP.h>
#include <gmp_lib/utils.h>
#include <cmath>

namespace as64_
{

namespace gmp_
{

class GMP_Update
{

public:

  typedef std::shared_ptr<GMP_Update> Ptr;

  /**GMP constructor.
   *  @param[in] gmp: n_DoF dmp.
   */
  GMP_Update(gmp_::GMP *gmp);

  void syncUpdate(bool set);

  void recursiveUpdate(bool set);

  void updateNow();

  void initSigmaw();

  void initSigmaWfromPosMsr(const arma::rowvec &x_data, double tol=1e-2);

  void initSigmaWfromVelMsr(const arma::rowvec &x_data, double Tf, double tol=1e-2);

  void initSigmaWfromAccelMsr(const arma::rowvec &x_data, double Tf, double tol=1e-2);

  void initExpSigmaw(double decay_rate=0.01);

  void setSigmaW(const arma::mat &Sw);

  arma::mat getSigmaW() const;

  void setMsrNoiseVar(double rv);

  // ==================================================
  // ===============   Online update  =================
  // ==================================================

  void updatePos(double x, const arma::vec &y, double r_n=nan(""));

  void updateVel(double x, double x_dot, const arma::vec &y_dot, double r_n=nan(""));

  void updateAccel(double x, double x_dot, double x_ddot, const arma::vec &y_ddot, double r_n=nan(""));


  /** Updates the weights based on 'n' measurements.
   *  @param[in] s: 1 x n vector of type gmp_::Phase, where the j-th entry is the phase for the j-th measurement.
   *  @param[in] Z: n_dof x n matrix, where the j-th column is the j-th measurement.
   *  @param[in] type: 1 x n vector of type gmp_::UPDATE_TYPE, where the j-th entry has the update type for the j-th measurement.
   *  @param[in] r_n: 1 x n vector, where the j-th entry is the noise variance for the j-th measurement (optional, default=this->rv).
   */
  void updateWeights(const std::vector<gmp_::Phase> &s, arma::mat Z, const std::vector<gmp_::UPDATE_TYPE> &type, arma::rowvec r_n={});

protected:

  gmp_::GMP *gmp; ///< pointer to nDoF GMP object

  arma::mat Sigma_w; ///< weights covariance for trajectory update
  double rv; ///< noise variance for trajectory update
  bool recursive_up;

  bool sync_up;

  void clearBatch();

  std::vector<gmp_::Phase> batch_s;
  arma::mat batch_Z;
  std::vector<gmp_::UPDATE_TYPE> batch_type;
  arma::rowvec batch_r_n;

}; // GMP_Update

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_UPDATE_H
