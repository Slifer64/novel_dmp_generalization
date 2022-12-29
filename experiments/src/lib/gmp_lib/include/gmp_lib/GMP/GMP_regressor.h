#ifndef GMP_LIB_GMP_REGRESSOR_H
#define GMP_LIB_GMP_REGRESSOR_H

// GMP_regressor class
// This class implements basic functionalities for a regressor vector
// constructed as a normalized sum of Gaussian kernel (basis) functions. 
// One can set the number of the kernels, their width and obtain the 
// regressor vector and its 1st, 2nd and 3rd derivatives.
// The kernels are equally spaced in [0 1] and take a one-dimensional
// input, say 'x', for example:
// phi = regressVec(x)
// phi_dot = regressVecDot(x, x_dot)
// phi_ddot = regressVecDDot(x, x_dot, x_ddot)
// phi_3dot = regressVec3Dot(x, x_dot, x_ddot, x_3dot)
//

#include <cmath>
#include <vector>
#include <cstring>
#include <memory>
#include <exception>
#include <cfloat>
#include <functional>
#include <fstream>
#include <armadillo>

namespace as64_
{

namespace gmp_
{

class GMP_regressor
{

// ====================================================
// ===============  Public functions  =================
// ====================================================

public:

  /** Constructor. Initialize with the number of kernels and their width.
  * @param[in] N_kernels The number of kernels (must be > 1).
  * @param[in] kernel_std_scaling Scaling of the kernel's std. (optional, default=1.0)
  */
  GMP_regressor(unsigned N_kernels, double kernel_std_scaling=1.0);

  /**
   * @brief Copy constructor.
   * @param[in] obj: Copy the contents of obj to *this. 
   */
  GMP_regressor(const GMP_regressor &obj);

  const GMP_regressor &operator=(const GMP_regressor &obj);

  /** Enable kernels truncation. 
  *  @param[in] zero_tol threshold below which the activation of a kernel is set to zero (optional, default=1e-8).
  *                      Set to zero, to disable kernels truncation.
  */
  void setTruncatedKernels(double zero_tol=1e-8);

  // ============================================================

  /** Returns the regressor vector.
  * @param[in] x The phase variable.
  * @return regressor vector.
  */
  arma::vec regressVec(double x) const;

  /** Returns the regressor vector 1st time derivative.
  * @param[in] x The phase variable.
  * @param[in] x_dot The phase variable 1st time derivative.
  * @return regressor vector 1st time derivative.
  */
  arma::vec regressVecDot(double x, double x_dot) const;

  /** Returns the regressor vector 2nd time derivative.
  * @param[in] x The phase variable.
  * @param[in] x_dot The phase variable 1st time derivative.
  * @param[in] x_ddot The phase variable 2nd time derivative.
  * @return regressor vector 2nd time derivative.
  */
  arma::vec regressVecDDot(double x, double x_dot, double x_ddot) const;

  /** Returns the regressor vector 3rd time derivative.
  * @param[in] x The phase variable.
  * @param[in] x_dot The phase variable 1st time derivative.
  * @param[in] x_ddot The phase variable 2nd time derivative.
  * @param[in] x_3dot The phase variable 3rd time derivative.
  * @return regressor vector 3rd time derivative.
  */
  arma::vec regressVec3Dot(double x, double dx, double ddx, double d3x) const;

  /** Returns the scaling of the kernels std.
  *  @return: the scaling of the kernels std.
  */
  double getKernelsStdScaling() const
  {
    return this->kernel_std_scaling;
  }

  /** Returns the centers of the kernels.
  *  @return: vector with the centers of the kernels.
  */
  arma::vec getCenters() const
  {
    return this->c;
  }

  /** Returns the inverse widths of the kernels.
  *  @return: vector with the invere width of the kernels.
  */
  arma::vec getInvWidths() const
  {
    return this->h;
  }

// =======================================================
// ===============  Protected functions  =================
// =======================================================
protected:

  // =============================================================

  /** Set the kernels centers and widths.
  * @param[in] N_kernels The number of kernels.
  * @param[in] kernel_std_scaling Scaling of the kernel's std.
  */
  void setKernels(unsigned N_kernels, double kernel_std_scaling);

  /** Set the kernels centers and widths.
  * @param[in] c Column vector with the center of each kernel.
  * @param[in] h Column vector with the width of each kernel.
  * @param[in] zero_tol Values equal or below this value are treated as zero (optional, default=DBL_MIN).
  */
  void setKernels(const arma::vec &c, const arma::vec &h, double zero_tol = DBL_MIN);

  /** Returns a column vector with the values of the kernel functions.
   * @param[in] x The phase variable.
   * @return: Column vector with the values of the kernel functions.
   */
  arma::vec kernelFun(double x) const;

  arma::vec truncGaussKernel(double x, double zero_tol);

  arma::vec kernelFunDot(double x, double dx) const;

  arma::vec kernelFunDDot(double x, double dx, double ddx) const;

  arma::vec kernelFun3Dot(double x, double dx, double ddx, double d3x) const;


// ========================================================
// ===============  Protected properties  =================
// ========================================================
protected:

  arma::vec c; ///< N_kernels x 1 vector with the kernels' centers
  arma::vec h; ///< N_kernels x 1 vector with the kernels' inverse width  

  std::function<arma::vec(double)> kernel_fun_ptr;

// ======================================================
// ===============  Private properties  =================
// ======================================================
private:

  // range of input x outside which psi = 0 due to finite numerical precision
  double x_min;
  double x_max; 
  
  double kernel_std_scaling; // the scaling of the kernel's std
  arma::vec kernel_std; // N_kernels x 1 vector with the kernels' std

  double zero_tol; // zero tolerance for kernels' truncation
};

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_REGRESSOR_H
