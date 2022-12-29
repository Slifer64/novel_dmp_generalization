#ifndef GMP_LIB_CONVERSIONS__64_H
#define GMP_LIB_CONVERSIONS__64_H


#include <armadillo>
#include <Eigen/Dense>

namespace as64_
{

namespace gmp_
{

inline void copyArmaToEigen(const arma::vec &a, Eigen::VectorXd &e)
{
  int n = a.size();
  e.resize(n);
  for (int i=0; i<n; i++) e(i) = a(i);
}

inline arma::vec getArmaFromEigen(const Eigen::VectorXd &e)
{
  int n = e.size();
  arma::vec a(n);
  for (int i=0; i<n; i++) a(i) = e(i);
  return a;
}


} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_CONVERSIONS__64_H