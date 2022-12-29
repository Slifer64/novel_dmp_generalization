#ifndef AS64_gmp_MATH_LIB_H
#define AS64_gmp_MATH_LIB_H

#include <Eigen/Dense>
#include <armadillo>

namespace as64_
{

namespace gmp_
{

Eigen::Vector4d quatInv(const Eigen::Vector4d &quat);
arma::vec quatInv(const arma::vec &quat);

Eigen::Matrix4d quat2mat(const Eigen::Vector4d &quat);

Eigen::Vector4d quatProd(const Eigen::Vector4d &quat1, const Eigen::Vector4d &quat2);
arma::vec quatProd(const arma::vec &quat1, const arma::vec &quat2);

Eigen::Vector4d quatExp(const Eigen::Vector3d &v_rot, double zero_tol=1e-16);
arma::vec quatExp(const arma::vec &v_rot, double zero_tol=1e-16);

Eigen::Vector3d quatLog(const Eigen::Vector4d &quat, double zero_tol=1e-16);
arma::vec quatLog(const arma::vec &quat, double zero_tol=1e-16);

arma::vec quatDiff(const arma::vec &quat1, const arma::vec &quat2);

arma::vec quatLogDiff(const arma::vec &quat1, const arma::vec &quat2);

// ===========================================================
// ===========================================================

/* Returns derivative of log given the rotational velocity and orientation (expressed w.r.t. the initial orientation)
 * @param[in] rotVel: Rotational velocity.
 * @param[in] Q1: Orientation expressed w.r.t. the initial orientation.
 * @return: Derivative of log.
 */
arma::vec rotVel_to_qLogDot(const arma::vec &rotVel, const arma::vec &logQ);

/** \brief TODO doc.*/
arma::vec qLogDot_to_rotVel(const arma::vec &logQ_dot, const arma::vec &Q);

/** \brief TODO doc.*/
arma::vec rotAccel_to_qLogDDot(const arma::vec &rotAccel, const arma::vec &rotVel, const arma::vec &Q);

/** \brief TODO doc.*/
arma::vec qLogDDot_to_rotAccel(const arma::vec &logQ_ddot, const arma::vec &rotVel, const arma::vec &Q);

/** \brief TODO doc.*/
arma::mat jacob_Q_qLog(const arma::vec &Q);

/** \brief TODO doc.*/
arma::mat jacob_qLog_Q(const arma::vec &Q);

/** \brief TODO doc.*/
arma::mat jacobDot_qLog_Q(const arma::vec &Q, const arma::vec &rotVel);

/** \brief TODO doc.*/
arma::mat jacobDot_Q_qLog(const arma::vec &Q, const arma::vec &rotVel);

} // namespace gmp_

} // namespace as64_

#endif // AS64_gmp_MATH_LIB_H
