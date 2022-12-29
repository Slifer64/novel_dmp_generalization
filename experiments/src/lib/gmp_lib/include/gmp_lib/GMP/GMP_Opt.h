#ifndef GMP_LIB_GMP_OPT_H
#define GMP_LIB_GMP_OPT_H

// N-DoF GMP Optimization class
//

#include <gmp_lib/GMP/GMP.h>

namespace as64_
{

namespace gmp_
{

class GMP_Opt
{
public:

  struct
  {
    double w_p; ///< position objective weight
    double w_v; ///< velocity objective weight
    double w_a; ///< acceleration objective weight

    long max_iter;
    bool warm_start;
    double time_limit;
    double abs_tol;
    double rel_tol;
  } settings;

  /** GMP constructor.
   * @param[in] gmp: n_DoF dmp.
   */
  GMP_Opt(gmp_::GMP *gmp);

  /** Trajectory optimization.
   * @param[in] num_points: Number of discreet points to use in the objective function.
   * @param[in] tau: Duration of motion.
   * @param[in] pos_constr: Vector of @GMPConstr position constraints. For no constraints pass '[]'.
   * @param[in] vel_constr: Vector of @GMPConstr velocity constraints. For no constraints pass '[]'.
   * @param[in] accel_constr: Vector of @GMPConstr acceleration constraints. For no constraints pass '[]'.
   * @return: 1 on success, 0 if suboptimal, -1 on failure.
   */
  int optimize(unsigned num_points = 200);

  int optimize(const arma::rowvec &x_data);

  void setMotionDuration(double tau);

  void setPosConstr(const arma::rowvec &x, const arma::mat &lb, const arma::mat &ub, const arma::rowvec &x_eq, const arma::mat &p_eq);

  void setVelConstr(const arma::rowvec &x, const arma::mat &lb, const arma::mat &ub, const arma::rowvec &x_eq, const arma::mat &v_eq);

  void setAccelConstr(const arma::rowvec &x, const arma::mat &lb, const arma::mat &ub, const arma::rowvec &x_eq, const arma::mat &a_eq);

  void clearPosConstr();

  void clearVelConstr();

  void clearAccelConstr();

  void clearConstr();

  std::string getExitMsg() const;

private:

  //ex_flag_map; ///< maps the exit flag value to the msg describing the exit flag
  std::string exit_msg;

  gmp_::GMP *gmp; ///< n_DoF GMP pointer

  double tau; ///< motion duration

  // Constraints
  arma::mat A_p;
  arma::mat pos_lb;
  arma::mat pos_ub;

  arma::mat Aeq_p;
  arma::mat pos_eq;

  arma::mat A_v;
  arma::mat vel_lb;
  arma::mat vel_ub;

  arma::mat Aeq_v;
  arma::mat vel_eq;

  arma::mat A_a;
  arma::mat accel_lb;
  arma::mat accel_ub;

  arma::mat Aeq_a;
  arma::mat accel_eq;

}; // GMP_Opt

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_OPT_H
