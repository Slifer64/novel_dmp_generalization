#ifndef GMP_LIB_MATH__64_H
#define GMP_LIB_MATH__64_H

#include <Eigen/Dense>
#include <armadillo>

namespace as64_
{

namespace gmp_
{

Eigen::Vector4d rotm2quat(Eigen::Matrix3d rotm, bool is_rotm_orthonormal=true);
arma::vec rotm2quat(const arma::mat &rotm, bool is_rotm_orthonormal=true);
Eigen::Vector4d mat2quat(Eigen::Matrix3d rotm);

Eigen::Matrix3d quat2rotm(Eigen::Vector4d quat);
arma::mat quat2rotm(const arma::vec &quat);

Eigen::Vector4d rotm2axang(Eigen::Matrix3d rotm);
arma::vec rotm2axang(const arma::mat &rotm);

Eigen::Matrix3d axang2rotm(Eigen::Vector4d axang);
arma::mat axang2rotm(const arma::vec &axang);

Eigen::Vector4d axang2quat(Eigen::Vector4d axang);
arma::vec axang2quat(const arma::vec &axang);

Eigen::Vector4d quat2axang(Eigen::Vector4d quat);
arma::vec quat2axang(const arma::vec &quat);

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_MATH__64_H
