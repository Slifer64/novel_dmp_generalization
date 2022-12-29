#ifndef GMP_LIB_GMP_MPC_H
#define GMP_LIB_GMP_MPC_H

// N-DoF GMP-MPC Optimization class
//

#include <gmp_lib/GMP/GMP.h>
#include <gmp_lib/math/linear_algebra.h>
#include <list>
#include <map>
#include <array>
#include <functional>

#include <armadillo>

#define GMP_MPC_L0G_STATISTICS

namespace as64_
{

namespace gmp_
{

class Obstacle
{
public:
  Obstacle(const std::string &name_): name(name_) {}
  virtual bool getConstrPlain(const arma::vec &p, const arma::vec &phi, arma::rowvec *Ai, double *bi) = 0;

  std::string getName() const { return name; }

public:
  arma::vec n_e;
  arma::vec p_e;
  arma::vec p;

private:
  std::string name;

};

class EllipseObstacle: public Obstacle
{
public:
  EllipseObstacle(const arma::vec &c_, const arma::mat &Sigma, const std::string &name_):
  Obstacle(name_), c(c_), inv_Sigma(arma::pinv(Sigma)) {}

  bool getConstrPlain(const arma::vec &p, const arma::vec &phi, arma::rowvec *Ai, double *bi) override
  {
    arma::vec &c = this->c;
    arma::mat inv_Sigma = 0.9*this->inv_Sigma; // enlarge to avoid numerical precision issues related to the solver's settings
    double temp = arma::as_scalar((p-c).t()*inv_Sigma*(p-c));

    if (temp > 1.1) return false;

    // std::cout << "p = " << p.t() << "\n";
    // std::cout << "c = " << c.t() << "\n";
    // std::cerr << "temp = " << temp << "\n";
    // std::cerr << "Sigma = \n" << arma::inv(inv_Sigma) << "\n";
    // exit(-1);

    // Find the point on the ellipsoid surface
    arma::vec p_e = c + (p-c) / std::sqrt(temp);
    // Find the norm to the ellipsoid on that point
    arma::vec n_e = inv_Sigma*(p_e - c);
    n_e /= arma::norm(n_e);
    arma::mat a_e = phi*n_e.t();
    // scale to avoid numeral issues related to the solver settings
    *Ai = 0.1*arma::vectorise(a_e).t();
    *bi = 0.1*arma::dot(n_e, p_e);

//    std::cout << "center: " << c.t() << "\n";
//    std::cout << "inv_Sigma: \n" << inv_Sigma << "\n";
//    std::cout << "Ai = \n";
//    for (int i=0; i<Ai->size(); i++) printf("%.4f  ", Ai->at(i));
//    printf("\nbi = %.3f\n", *bi);
//    exit(-1);

    this->p_e = p_e;
    this->n_e = n_e;
    this->p = p;

    return true;
  }

private:
  arma::vec c;
  arma::mat inv_Sigma;
};

class GMP_MPC
{

public:

  typedef std::shared_ptr<GMP_MPC> Ptr;

  struct Log
  {
    void clear()
    { n_e_data.clear(); p_e_data.clear(); p_data.clear(); si_data.clear(); yd_points.clear(); y_pred_points.clear(); }

    std::vector<arma::vec> n_e_data;
    std::vector<arma::vec> p_e_data;
    std::vector<arma::vec> p_data;
    std::vector<double> si_data;
    std::vector<arma::vec> yd_points;
    std::vector<arma::vec> y_pred_points;
    arma::vec y_current;
    arma::vec yd_current;
    arma::vec y_target;
    arma::mat W_opt;
    arma::mat Wd;
  } log;

  std::function<void(const Log &)> plot_callback;

  struct
  {
    unsigned max_iter;
    double time_limit;
    double abs_tol;
    double rel_tol;
  } settings;

  struct Solution
  {
    arma::vec y;
    arma::vec y_dot;
    arma::vec y_ddot;

    arma::vec pos_slack;
    arma::vec vel_slack;
    arma::vec accel_slack;

    int exit_flag;
    std::string exit_msg;
  };
  
  GMP_MPC(const GMP *gmp, unsigned N_horizon, double pred_time_step, unsigned N_kernels, 
          double kernel_std_scaling, const std::array<double,3> &slack_gains, double trunc_kern_thres=1e-6);

  void setInitialState(const arma::vec &y0, const arma::vec &y0_dot, const arma::vec &y0_ddot, double s, double s_dot, double s_ddot);
  
  void setFinalState(const arma::vec &yf, const arma::vec &yf_dot, const arma::vec &yf_ddot, double s, double s_dot, double s_ddot);
  void setFinalState(const arma::vec &yf, const arma::vec &yf_dot, const arma::vec &yf_ddot, double s, double s_dot, double s_ddot, const arma::vec &err_tol);

  void setObjCostGains(double pos_gain, double vel_gain);

  void setPosLimits(const arma::vec &lb, const arma::vec &ub);
  void setPosLimits(const Eigen::VectorXd &lb, const Eigen::VectorXd &ub);
  
  void setVelLimits(const arma::vec &lb, const arma::vec &ub);
  void setVelLimits(const Eigen::VectorXd &lb, const Eigen::VectorXd &ub);
  
  void setAccelLimits(const arma::vec &lb, const arma::vec &ub);
  void setAccelLimits(const Eigen::VectorXd &lb, const Eigen::VectorXd &ub);
  
  void addViaPoint(double s, const arma::vec &y, double err_tol = 1e-6);

  void addEllipsoidObstacle(const arma::vec &c, const arma::mat &Sigma, std::string name="");
  void removeObstacle(const std::string &name);

  void setPosSlackLimit(double s_lim);
  
  void setVelSlackLimit(double s_lim);
  
  void setAccelSlackLimit(double s_lim);

  void setCanonicalSystemFunction(std::function<std::array<double,2>(double, double)> can_sys_fun);

  GMP_MPC::Solution solve(double s, double s_dot);

  arma::vec getYd(double s) const;

  arma::vec getYdDot(double s, double s_dot) const;

  arma::vec getYdDDot(double s, double s_dot, double s_ddot) const;

  std::vector<double> elaps_t_data;
  std::vector<double> calcH_elaps_t_data;
  std::vector<double> solve_elaps_t_data;

  #ifdef GMP_MPC_L0G_STATISTICS
    Eigen::SpMat H_;
    Eigen::SpMat A_;
    double s_rec = 0.5;
  #endif

protected:

  void process_obstacles(const arma::vec &p1, const arma::vec &phi, arma::rowvec *Ai, double *bi);

  struct Ellipsoid
  {
    Ellipsoid() {}
    Ellipsoid(const arma::vec &c_, const arma::vec &Sigma_): c(c_), inv_Sigma(arma::pinv(Sigma_)) {}
    arma::vec c;
    arma::mat inv_Sigma;
  };
//  std::map<std::string, Ellipsoid> elips_obst;
  std::map<std::string, std::shared_ptr<Obstacle>> obst_map;
  unsigned n_obst;

  void integratePhase(double s, double s_dot, std::vector<double> &si_data, std::vector<double> &si_dot_data, std::vector<double> &si_ddot_data);

  Eigen::SpMat getAjMat(double s, double s_dot, double s_ddot) const;

  #ifdef GMP_MPC_L0G_STATISTICS
    arma::wall_clock timer;
    arma::wall_clock timer_H;
    arma::wall_clock timer_solve;
  #endif

  std::function<std::array<double,2>(double, double)> can_sys_fun;

  struct ViaPoint
  {
    ViaPoint(double s_, const Eigen::VectorXd &pos_, double err_tol_):
    s(s_), pos(pos_), err_tol(err_tol_) {}

    double s;
    Eigen::VectorXd pos;
    Eigen::SpMat Phi;
    double err_tol;

    Eigen::VectorXd z_dual;
  };
  std::list<ViaPoint> via_points;
  unsigned n_via; // for how many via-points to allocate space
  
  arma::mat W_mpc; ///< mpc optimized weights
  const GMP *gmp_ref;

  GMP::Ptr gmp_mpc;

  unsigned N_kernels;

  unsigned N;
  std::vector<double> dt_;

  unsigned n_ineq; // number of inequality constraints
  unsigned n_var; // number of optimization variables
  unsigned n_eq; // number of equality constraints
  
  Eigen::VectorXd Qs_vec;
  std::vector<Eigen::Triplet<double>> Aineq_slack_triplets;
  
  double sqrt_pos_gain;
  double sqrt_vel_gain;
  Eigen::VectorXd sqrt_pos_vel_gain_vec;

  // Eigen::SpMat Qi;
  // Eigen::SpMat QN;

  arma::vec Z0;
  arma::vec Z0_dual_ineq;
  arma::vec Z0_dual_eq;
  
  Eigen::VectorXd pos_lb;
  Eigen::VectorXd pos_ub;
  
  Eigen::VectorXd vel_lb;
  Eigen::VectorXd vel_ub;
  
  Eigen::VectorXd accel_lb;
  Eigen::VectorXd accel_ub;
  
  unsigned n_slack;
  
  bool pos_slack;
  bool vel_slack;
  bool accel_slack;
  
  Eigen::VectorXd pos_slack_lim;
  Eigen::VectorXd vel_slack_lim;
  Eigen::VectorXd accel_slack_lim;
  
  Eigen::SpMat Phi0;
  Eigen::VectorXd x0;
  
  double s_f;
  std::vector<Eigen::Triplet<double>> Phi_f_triplets;
  Eigen::SpMat Phi_f;
  Eigen::SpMat Qf;
  arma::vec phi_f, phi_f_dot, phi_f_ddot;
  Eigen::VectorXd x_f;
  Eigen::VectorXd err_tol_f;

private:

  unsigned n_dof;
  unsigned n_dof3;

  double inf;

  Eigen::SpMat I_ndof;
  Eigen::VectorXd inf_ndof;
  Eigen::VectorXd ones_ndof;
  Eigen::VectorXd zeros_ndof;

    
}; // class GMP_MPC

} // namespace gmp_

} // namespace as64_

#endif // GMP_LIB_GMP_MPC_H
