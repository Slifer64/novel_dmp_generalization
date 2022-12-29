#include <gmp_lib/GMP/GMP_regressor.h>

namespace as64_
{

namespace gmp_
{

// public:

  GMP_regressor::GMP_regressor(unsigned N_kernels, double kernel_std_scaling)
  {
    this->setTruncatedKernels(0);
    this->setKernels(N_kernels, kernel_std_scaling);
  }

  void GMP_regressor::setTruncatedKernels(double zero_tol)
  {
    if (zero_tol < 0) throw std::runtime_error("Kernels truncation threshold must be non-negative.");

    this->zero_tol = zero_tol;
    if (zero_tol)
    {
      setKernels(this->c, this->h, zero_tol);
      this->kernel_fun_ptr = std::bind(&GMP_regressor::truncGaussKernel,this, std::placeholders::_1, zero_tol);
    }
    else
      this->kernel_fun_ptr = std::bind(&GMP_regressor::kernelFun,this,std::placeholders::_1);
    
  }

  GMP_regressor::GMP_regressor(const GMP_regressor &obj)
  {
    *this = obj;
  }

  const GMP_regressor &GMP_regressor::operator=(const GMP_regressor &obj)
  {
    this->c = obj.c;
    this->h = obj.h;
    this->x_min = obj.x_min;
    this->x_max = obj.x_max;
    this->kernel_std_scaling = obj.kernel_std_scaling;
    this->kernel_std = obj.kernel_std;
    this->zero_tol = obj.zero_tol;
    this->setTruncatedKernels(obj.zero_tol);

    return *this;
  }

  // ============================================================

  arma::vec GMP_regressor::regressVec(double x) const
  {
    arma::vec psi(this->c.size());

    // take appropriate actions when x causes phi = 0 due to finite
    // numerical precision.
    if (x < this->x_min)
    {
      psi = arma::vec().zeros(psi.size());
      psi(0) = 1;
    }
    else if (x > this->x_max)
    {
      psi = arma::vec().zeros(psi.size());
      psi.back() = 1;
    }
    else
    {
      psi = this->kernel_fun_ptr(x);
    }

    // if ( arma::sum(psi)  < DBL_MIN)
    // {
    //   std::cerr << "x=" << x << " : ZERO! << (" << this->x_min << ", " << this->x_max << ")\n";
    // }
    
    return psi / arma::sum(psi);
  }

  arma::vec GMP_regressor::regressVecDot(double x, double x_dot) const
  {
    if (x < this->x_min || x > this->x_max) return arma::vec().zeros(this->c.size());

    arma::vec psi = this->kernel_fun_ptr(x);
    arma::vec psi_dot = this->kernelFunDot(x, x_dot);
    double sum_psi = arma::sum(psi);
    double sum_psi_dot = arma::sum(psi_dot);

    arma::vec phi = psi / sum_psi;
    return ( psi_dot - phi*sum_psi_dot ) / sum_psi;
  }

  arma::vec GMP_regressor::regressVecDDot(double x, double x_dot, double x_ddot) const
  {
    if (x < this->x_min || x > this->x_max) return arma::vec().zeros(this->c.size());

    arma::vec psi = this->kernel_fun_ptr(x);
    arma::vec psi_dot = this->kernelFunDot(x, x_dot);
    arma::vec psi_ddot = this->kernelFunDDot(x, x_dot, x_ddot);
    double sum_psi = arma::sum(psi);
    double sum_psi_dot = arma::sum(psi_dot);
    double sum_psi_ddot = arma::sum(psi_ddot);

    arma::vec phi = psi / sum_psi;
    arma::vec phi_dot = ( psi_dot - phi*sum_psi_dot ) / sum_psi;
    return (psi_ddot - 2*phi_dot*sum_psi_dot - phi*sum_psi_ddot) / sum_psi;
  }

  arma::vec GMP_regressor::regressVec3Dot(double x, double dx, double ddx, double d3x) const
  {
    if (x < this->x_min || x > this->x_max) return arma::vec().zeros(this->c.size());

    arma::vec psi = this->kernel_fun_ptr(x);
    arma::vec psi_dot = this->kernelFunDot(x, dx);
    arma::vec psi_ddot = this->kernelFunDDot(x, dx, ddx);
    arma::vec psi_3dot = this->kernelFun3Dot(x, dx, ddx, d3x);
    double sum_psi = arma::sum(psi);
    double sum_psi_dot = arma::sum(psi_dot);
    double sum_psi_ddot = arma::sum(psi_ddot);
    double sum_psi_3dot = arma::sum(psi_3dot);

    arma::vec phi = psi / sum_psi;
    arma::vec phi_dot = ( psi_dot - phi*sum_psi_dot ) / sum_psi;
    arma::vec phi_ddot = (psi_ddot - 2*phi_dot*sum_psi_dot - phi*sum_psi_ddot) / sum_psi;
    return (psi_3dot - 3*phi_ddot*sum_psi_dot - 3*phi_dot*sum_psi_ddot - phi*sum_psi_3dot) / sum_psi;
  }

// protected:

  // =============================================================

  arma::vec GMP_regressor::kernelFun(double x) const
  {
    return arma::exp(-this->h % (arma::pow(x-this->c,2)));
  }

  arma::vec GMP_regressor::truncGaussKernel(double x, double zero_tol)
  {
    arma::vec psi = this->kernelFun(x);
    arma::uvec ind = arma::find(psi < zero_tol);
    psi.elem(ind) = arma::vec().zeros(ind.size());
    return psi;
  }

  arma::vec GMP_regressor::kernelFunDot(double x, double dx) const
  {
    arma::vec psi = this->kernel_fun_ptr(x);
    arma::vec a = (x-this->c)*dx;
    return -2*this->h%( psi % a );
  }

  arma::vec GMP_regressor::kernelFunDDot(double x, double dx, double ddx) const
  {
    arma::vec psi = this->kernel_fun_ptr(x);
    arma::vec psi_dot = this->kernelFunDot(x, dx);

    arma::vec a = (x-this->c)*dx;
    arma::vec a_dot = (x-this->c)*ddx + std::pow(dx,2);
    return -2*this->h%( psi_dot%a + psi%a_dot );
  }

  arma::vec GMP_regressor::kernelFun3Dot(double x, double dx, double ddx, double d3x) const
  {
    arma::vec psi = this->kernel_fun_ptr(x);
    arma::vec psi_dot = this->kernelFunDot(x, dx);
    arma::vec psi_ddot = this->kernelFunDDot(x, dx, ddx);

    arma::vec a = (x-this->c)*dx;
    arma::vec a_dot = (x-this->c)*ddx + std::pow(dx,2);
    arma::vec a_ddot = (x-this->c)*d3x + 3*dx*ddx;
    return -2*this->h%( psi_ddot%a + 2*psi_dot%a_dot + psi%a_ddot );
  }

  void GMP_regressor::setKernels(unsigned N_kernels, double kernel_std_scaling)
  {
    if (N_kernels < 2) throw std::runtime_error("[GMP_regressor::setKernels]: At least 2 kernels are required!");

    arma::vec c = arma::linspace<arma::vec>(0,1, N_kernels);

    double hi = 1 / std::pow(kernel_std_scaling*(c(1)-c(0)),2);
    arma::vec h = arma::vec().ones(N_kernels) * hi;

    setKernels(c, h);
  }

  void GMP_regressor::setKernels(const arma::vec &c, const arma::vec &h, double zero_tol)
  {
    if (c.size() < 2) throw std::runtime_error("[GMP_regressor::setKernels]: At least 2 kernels are required!");

    this->c = c;
    this->h = h;

    // assuming all kernels have the same width
    this->kernel_std_scaling = std::sqrt(1.0 / h(0)) / ( c(1) - c(0) );

    // find range of x outside which psi = 0 due to finite numerical precision
    this->x_min = this->c(0) - std::sqrt( -std::log(zero_tol) / this->h(0) );
    this->x_max = this->c.back() + std::sqrt( -std::log(zero_tol) / this->h.back() );
  }

} // namespace gmp_

} // namespace as64_
