#include <gmp_lib/GMP/GMP_MPC.h>
#include <gmp_lib/math/conversions.h>

#include <cfloat> // for infinity

//#include <osqp_lib/osqp.h>
//#include <osqp_lib/constants.h>
//#include <osqp_lib/csc_mat.h>
#include <osqp_lib/qp_solver.h>

#include <boost/numeric/odeint.hpp>
// using namespace boost::numeric::odeint;


#include <iomanip>

namespace as64_
{

namespace gmp_
{

#define GMP_MPC_fun_ std::string("[GMP_MPC::") + __func__ + "]: "

typedef std::array<double,2> ode_state;
// typedef boost::numeric::odeint::euler< ode_state > euler_stepper_type;
typedef boost::numeric::odeint::runge_kutta_cash_karp54< ode_state > error_stepper_type;
typedef boost::numeric::odeint::controlled_runge_kutta< error_stepper_type > controlled_stepper_type;

arma::mat blkdiag(const arma::mat &A, const arma::mat &B)
{
  if (A.empty()) return B;

  arma::mat C(A.n_rows + B.n_rows, A.n_cols + B.n_cols);
  C.submat(0,0, A.n_rows-1, A.n_cols-1) = A;
  C.submat(A.n_rows, A.n_cols, C.n_rows-1, C.n_cols-1) = B;
  return C;
}

arma::mat blkdiag(const arma::mat &A, const arma::mat &B, const arma::mat &C)
{
  return blkdiag( blkdiag(A, B), C );
}

arma::mat blkdiag(const arma::mat &A, double b)
{
  arma::mat C(A.n_rows + 1, A.n_cols + 1);
  if (!A.empty()) C.submat(0,0, A.n_rows-1, A.n_cols-1) = A;
  C.submat(A.n_rows, A.n_cols, C.n_rows-1, C.n_cols-1) = b;
  return C;
}

void addTriplets(std::vector<Eigen::Triplet<double>> &a, const std::vector<Eigen::Triplet<double>> &b,
                 unsigned row_offset=0, unsigned col_offset=0)
{
  for (const Eigen::Triplet<double> &t : b)
    a.push_back( Eigen::Triplet<double>( t.row()+row_offset , t.col()+col_offset, t.value() ) );
}

void addTriplets(std::vector<Eigen::Triplet<double>> &a, const arma::mat &m,
                 unsigned row_offset=0, unsigned col_offset=0)
{
  for (int i=0; i<m.n_rows; i++)
  {
    for (int j=0; j<m.n_cols; j++)
    {
      if (m(i,j)) a.push_back( Eigen::Triplet<double>( i+row_offset , j+col_offset, m(i, j) ) );
    }
  }
}

void addDiagonalTriplets(std::vector<Eigen::Triplet<double>> &coeff, const Eigen::VectorXd &v,
                         unsigned row_offset=0, unsigned col_offset=0)
{
  for (int i=0; i<v.size(); i++) coeff.push_back(Eigen::Triplet<double>(row_offset+i,col_offset+i,v[i]));
}

void addDiagonalTriplets(std::vector<Eigen::Triplet<double>> &coeff, unsigned dim, double value,
                         unsigned row_offset=0, unsigned col_offset=0)
{
  for (int i=0; i<dim; i++) coeff.push_back(Eigen::Triplet<double>(row_offset+i,col_offset+i,value));
}

void addToAineqTriplets(std::vector<Eigen::Triplet<double>> &Aineq_triplets, unsigned n_dof,
                        const arma::vec &phi, const arma::vec &phi_dot, const arma::vec &phi_ddot,
                        unsigned row_offset)
{
  int N_kernels = phi.size();

  for (int j=0; j<N_kernels; j++)
  {
    if (phi(j))
    {
      for (int i=0; i<n_dof; i++) Aineq_triplets.push_back( Eigen::Triplet<double>( row_offset+i, j+i*N_kernels, phi(j) ) );
    }
    if (phi_dot(j))
    {
      for (int i=0; i<n_dof; i++) Aineq_triplets.push_back( Eigen::Triplet<double>( row_offset+n_dof+i, j+i*N_kernels, phi_dot(j) ) );
    }
    if (phi_ddot(j))
    {
      for (int i=0; i<n_dof; i++) Aineq_triplets.push_back( Eigen::Triplet<double>( row_offset+2*n_dof+i, j+i*N_kernels, phi_ddot(j) ) );
    }
  }

}

Eigen::SpMat createPsiMat(int n_dof, const arma::vec &phi, const arma::vec &phi_dot)
{
  int N_kernels = phi.size();

  std::vector<Eigen::Triplet<double>> coeff;
  coeff.reserve(N_kernels*n_dof*2);  // upper bound of nonzero elements
  for (int j=0; j<N_kernels; j++)
  {
    if (phi(j))
    {
      for (int i=0; i<n_dof; i++) coeff.push_back( Eigen::Triplet<double>( i, j+i*N_kernels, phi(j) ) );
    }
    if (phi_dot(j))
    {
      for (int i=0; i<n_dof; i++) coeff.push_back( Eigen::Triplet<double>( n_dof+i, j+i*N_kernels, phi_dot(j) ) );
    }
  }

  Eigen::SpMat sm(3*n_dof, n_dof*N_kernels);
  sm.setFromTriplets(coeff.begin(), coeff.end());
  return sm;
}

Eigen::SpMat GMP_MPC::getAjMat(double s, double s_dot, double s_ddot) const
{
  arma::vec phi = this->gmp_mpc->regressVec(s);
  arma::vec phi_dot = this->gmp_mpc->regressVecDot(s, s_dot);
  arma::vec phi_ddot = this->gmp_mpc->regressVecDDot(s, s_dot, s_ddot);

  std::vector<Eigen::Triplet<double>> coeff;
  coeff.reserve( 3*n_dof*N_kernels ); // upper bound of nonzero elements
  for (int j=0; j<N_kernels; j++)
  {
    if (phi(j))
    {
      for (int i=0; i<n_dof; i++) coeff.push_back( Eigen::Triplet<double>( i, j+i*N_kernels, phi(j) ) );
    }
    if (phi_dot(j))
    {
      for (int i=0; i<n_dof; i++) coeff.push_back( Eigen::Triplet<double>( n_dof+i, j+i*N_kernels, phi_dot(j) ) );
    }
    if (phi_ddot(j))
    {
      for (int i=0; i<n_dof; i++) coeff.push_back( Eigen::Triplet<double>( 2*n_dof+i, j+i*N_kernels, phi_ddot(j) ) );
    }
  }

  Eigen::SpMat sm(3*n_dof, n_dof*N_kernels);
  sm.setFromTriplets(coeff.begin(), coeff.end());
  return sm;
}

// ====================================================================
// ====================================================================
// ====================================================================

GMP_MPC::GMP_MPC(const GMP *gmp, unsigned N_horizon, double pred_time_step, unsigned N_kernels, 
                double kernel_std_scaling, const std::array<double,3> &slack_gains, double trunc_kern_thres)
{
  this->inf = OSQP_INFTY;

  this->settings.time_limit = 0;
  this->settings.max_iter = 12000;
  this->settings.abs_tol = 1e-3;
  this->settings.rel_tol = 1e-5;

  this->n_dof = gmp->numOfDoFs();
  this->n_dof3 = 3*n_dof;
  this->I_ndof = Eigen::SpMat(n_dof,n_dof);
  this->I_ndof.setIdentity();
  
  this->inf_ndof = Eigen::VectorXd(n_dof).setConstant(this->inf);
  this->ones_ndof = Eigen::VectorXd(n_dof).setOnes();
  this->zeros_ndof = Eigen::VectorXd(n_dof).setZero();

  this->gmp_ref = gmp;
  
  this->gmp_mpc.reset(new GMP(n_dof, N_kernels, kernel_std_scaling));
  this->gmp_mpc->setTruncatedKernels(trunc_kern_thres);

  this->W_mpc = arma::mat().zeros(n_dof, N_kernels);

  this->N_kernels = N_kernels;
  
  this->N = N_horizon;
  this->dt_.resize(this->N);
  for (int i=0; i<dt_.size(); i++) this->dt_[i] = pred_time_step;

  this->pos_slack = slack_gains[0] > 0;
  this->vel_slack = slack_gains[1] > 0;
  this->accel_slack = slack_gains[2] > 0;
  this->n_slack = n_dof * (this->pos_slack + this->vel_slack + this->accel_slack);

  this->n_via = 3; // number of via-points to consider at each optimization window

  this->n_var = n_dof*N_kernels + n_slack;
  unsigned n_bound = N * n_dof3; // pos,vel,accel limits
  unsigned n_final_state = n_dof3;
  this->n_obst = N;
  unsigned n_via_constr = n_via*n_dof3;
  this->n_ineq = n_bound + n_final_state + n_obst + n_via_constr + n_slack;
  this->n_eq = n_dof3; // initial state constraint

  this->Z0 = arma::vec().zeros(n_var);
  this->Z0_dual_ineq = arma::vec().zeros(n_ineq);
  this->Z0_dual_eq = arma::vec().zeros(n_eq);

  Eigen::VectorXd Qp_s, Qv_s, Qa_s;
  // Eigen::SpMat Ap_s, Av_s, Aa_s;
  Aineq_slack_triplets.reserve(3*n_dof);
  if (this->pos_slack)
  {
    Qp_s = slack_gains[0]*ones_ndof;
    for (int i=0; i<n_dof; i++) Aineq_slack_triplets.push_back(Eigen::Triplet<double>(i,i, -1.0));
    //Ap_s = Eigen::join_vert( -I_ndof, Eigen::SpMat(2*n_dof, n_dof) ); 
  }
  if (this->vel_slack)
  {
    Qv_s = slack_gains[1]*ones_ndof;
    for (int i=n_dof; i<2*n_dof; i++) Aineq_slack_triplets.push_back(Eigen::Triplet<double>(i,i, -1.0));
    // Av_s = Eigen::join_vert( Eigen::SpMat(n_dof, n_dof), -I_ndof, Eigen::SpMat(n_dof, n_dof) );
  }
  if (this->accel_slack)
  {
    Qa_s = slack_gains[2]*ones_ndof;
    for (int i=2*n_dof; i<3*n_dof; i++) Aineq_slack_triplets.push_back(Eigen::Triplet<double>(i,i, -1.0));
    // Aa_s = Eigen::join_vert( Eigen::SpMat(2*n_dof, n_dof), -I_ndof );
  }

  this->Qs_vec.resize(Qp_s.size() + Qv_s.size() + Qa_s.size());
  this->Qs_vec << Qp_s, Qv_s, Qa_s;

  setObjCostGains(1.0, 0.0);

  this->setPosLimits(-inf_ndof, inf_ndof);
  this->setVelLimits(-inf_ndof, inf_ndof);
  this->setAccelLimits(-inf_ndof, inf_ndof);

  this->setPosSlackLimit(this->inf);
  this->setVelSlackLimit(this->inf);
  this->setAccelSlackLimit(this->inf);

  Qf.resize(n_dof3, n_dof3);
  Qf.setIdentity();
  Eigen::VectorXd vf(n_dof3);
  vf << 1e3*ones_ndof, 1e2*ones_ndof, 1e1*ones_ndof;
  Qf.diagonal() = vf;
}

void GMP_MPC::setObjCostGains(double pos_gain, double vel_gain)
{
  if (pos_gain < 0) throw std::runtime_error(GMP_MPC_fun_ + " 'pos_gain' must be non-negative.");
  if (vel_gain < 0) throw std::runtime_error(GMP_MPC_fun_ + " 'vel_gain' must be non-negative."); 

  this->sqrt_pos_gain = std::sqrt(pos_gain);
  this->sqrt_vel_gain = std::sqrt(vel_gain);

  sqrt_pos_vel_gain_vec.resize(2*n_dof);
  sqrt_pos_vel_gain_vec << Eigen::VectorXd::Constant(n_dof, sqrt_pos_gain), Eigen::VectorXd::Constant(n_dof, sqrt_vel_gain);
  
  // this->Qi = Eigen::SpMat(2*n_dof,2*n_dof);
  // this->Qi.setIdentity();
  // this->Qi.diagonal() = v;

  // // this->Qi = blkdiag( pos_gain*I_ndof, vel_gain*I_ndof );
  // this->QN = this->Qi;
}

void GMP_MPC::setInitialState(const arma::vec &y0, const arma::vec &y0_dot, const arma::vec &y0_ddot, double s, double s_dot, double s_ddot)
{
  if (y0.size() != this->n_dof) throw std::runtime_error(GMP_MPC_fun_ + "y0.size() must be equal to n_dof");
  if (y0_dot.size() != this->n_dof) throw std::runtime_error(GMP_MPC_fun_ + "y0_dot.size() must be equal to n_dof");
  if (y0_ddot.size() != this->n_dof) throw std::runtime_error(GMP_MPC_fun_ + "y0_ddot.size() must be equal to n_dof");

  arma::vec phi0_ = this->gmp_mpc->regressVec(s);
  arma::vec phi0_dot_ = this->gmp_mpc->regressVecDot(s, s_dot);
  arma::vec phi0_ddot_ = this->gmp_mpc->regressVecDDot(s, s_dot, s_ddot);
  
  this->Phi0 = getAjMat(s, s_dot, s_ddot);
  this->x0.resize(3*n_dof);
  this->x0 << Eigen::Map<const Eigen::VectorXd>(y0.memptr(), n_dof),
              Eigen::Map<const Eigen::VectorXd>(y0_dot.memptr(), n_dof),
              Eigen::Map<const Eigen::VectorXd>(y0_ddot.memptr(), n_dof);
}

void GMP_MPC::setFinalState(const arma::vec &yf, const arma::vec &yf_dot, const arma::vec &yf_ddot, double s, double s_dot, double s_ddot)
{
  setFinalState(yf, yf_dot, yf_ddot, s, s_dot, s_ddot, arma::vec().zeros(3) );
}

void GMP_MPC::setFinalState(const arma::vec &yf, const arma::vec &yf_dot, const arma::vec &yf_ddot, double s, double s_dot, double s_ddot, const arma::vec &err_tol)
{

  if (yf.size() != this->n_dof) throw std::runtime_error(GMP_MPC_fun_ + "yf.size() must be equal to n_dof");
  if (yf_dot.size() != this->n_dof) throw std::runtime_error(GMP_MPC_fun_ + "yf_dot.size() must be equal to n_dof");
  if (yf_ddot.size() != this->n_dof) throw std::runtime_error(GMP_MPC_fun_ + "yf_ddot.size() must be equal to n_dof");
  if (err_tol.size() != 3) throw std::runtime_error(GMP_MPC_fun_ + "err_tol.size() must be equal to 3");

  this->s_f = s;

  this->err_tol_f.resize(3*n_dof);
  this->err_tol_f << err_tol(0)*ones_ndof, err_tol(1)*ones_ndof, err_tol(2)*ones_ndof;

  this->phi_f = this->gmp_mpc->regressVec(s);
  this->phi_f_dot = this->gmp_mpc->regressVecDot(s, s_dot);
  this->phi_f_ddot = this->gmp_mpc->regressVecDDot(s, s_dot, s_ddot);

  Phi_f_triplets.clear();
  Phi_f_triplets.reserve(N_kernels*n_dof*3);
  for (int j=0; j<N_kernels; j++)
  {
    if (phi_f(j))
    {
      for (int i=0; i<n_dof; i++) Phi_f_triplets.push_back( Eigen::Triplet<double>( i, j+i*N_kernels, phi_f(j) ) );
    }
    if (phi_f_dot(j))
    {
      for (int i=0; i<n_dof; i++) Phi_f_triplets.push_back( Eigen::Triplet<double>( n_dof+i, j+i*N_kernels, phi_f_dot(j) ) );
    }
    if (phi_f_ddot(j))
    {
      for (int i=0; i<n_dof; i++) Phi_f_triplets.push_back( Eigen::Triplet<double>( 2*n_dof+i, j+i*N_kernels, phi_f_ddot(j) ) );
    }
  }

  Phi_f.resize(3*n_dof, N_kernels*n_dof);
  Phi_f.setFromTriplets(Phi_f_triplets.begin(), Phi_f_triplets.end());

  this->x_f.resize(3*n_dof);
  this->x_f << Eigen::Map<const Eigen::VectorXd>(yf.memptr(), n_dof),
               Eigen::Map<const Eigen::VectorXd>(yf_dot.memptr(), n_dof),
               Eigen::Map<const Eigen::VectorXd>(yf_ddot.memptr(), n_dof);
}

void GMP_MPC::addViaPoint(double s, const arma::vec &y, double err_tol)
{    
  if (y.size() != this->n_dof) throw std::runtime_error(GMP_MPC_fun_ + "y.size() must be equal to n_dof");

  ViaPoint vp(s, Eigen::Map<const Eigen::VectorXd>(y.memptr(),y.size()), err_tol);


  arma::vec phi = gmp_mpc->regressVec(s);
  std::vector<Eigen::Triplet<double>> coeff;
  coeff.reserve(N_kernels*n_dof);
  for (int j=0; j<N_kernels; j++)
  {
    if (phi(j))
    {
      for (int i=0; i<n_dof; i++) coeff.push_back( Eigen::Triplet<double>( i, j+i*N_kernels, phi(j) ) );
    }
  }

  vp.Phi.resize(n_dof, n_dof*N_kernels);
  vp.Phi.setFromTriplets(coeff.begin(), coeff.end());

  vp.z_dual = Eigen::VectorXd::Zero(this->n_dof3);

  // find where to insert the new via-point, so that they are sorted
  unsigned n = via_points.size();
  auto it = via_points.rbegin();
  while (it != via_points.rend()) // search backwards, since new via-points will probably refer to future time
  {
    if (s > it->s) break;
    it++;
  }   
  via_points.insert(it.base(), vp); // base returns one element after 'it'
}

void GMP_MPC::addEllipsoidObstacle(const arma::vec &c, const arma::mat &Sigma, std::string name)
{
  if (name.empty()) name = std::to_string(this->obst_map.size());
  this->obst_map[name] = std::shared_ptr<Obstacle>(new EllipseObstacle(c, Sigma, name));
}

void GMP_MPC::removeObstacle(const std::string &name)
{
  auto it = this->obst_map.find(name);
  if (it != this->obst_map.end()) this->obst_map.erase(it);
}

// -----------------------------------------------
    
void GMP_MPC::setPosLimits(const arma::vec &lb, const arma::vec &ub)
{
  Eigen::Map<const Eigen::VectorXd> lb_(lb.memptr(), lb.size());
  Eigen::Map<const Eigen::VectorXd> ub_(ub.memptr(), ub.size());
  setPosLimits(lb_, ub_);
}

void GMP_MPC::setPosLimits(const Eigen::VectorXd &lb, const Eigen::VectorXd &ub)
{
  if (lb.size() != this->n_dof) throw std::runtime_error(GMP_MPC_fun_ + "lb.size() must be equal to n_dof");
  if (ub.size() != this->n_dof) throw std::runtime_error(GMP_MPC_fun_ + "ub.size() must be equal to n_dof");

  this->pos_lb = lb;
  this->pos_ub = ub;
}

// -----------------------------------------------

void GMP_MPC::setVelLimits(const arma::vec &lb, const arma::vec &ub)
{
  Eigen::Map<const Eigen::VectorXd> lb_(lb.memptr(), lb.size());
  Eigen::Map<const Eigen::VectorXd> ub_(ub.memptr(), ub.size());
  setVelLimits(lb_, ub_);
}

void GMP_MPC::setVelLimits(const Eigen::VectorXd &lb, const Eigen::VectorXd &ub)
{
  if (lb.size() != this->n_dof) throw std::runtime_error(GMP_MPC_fun_ + "lb.size() must be equal to n_dof");
  if (ub.size() != this->n_dof) throw std::runtime_error(GMP_MPC_fun_ + "ub.size() must be equal to n_dof");

  this->vel_lb = lb;
  this->vel_ub = ub;
}

// -----------------------------------------------

void GMP_MPC::setAccelLimits(const arma::vec &lb, const arma::vec &ub)
{ 
  Eigen::Map<const Eigen::VectorXd> lb_(lb.memptr(), lb.size());
  Eigen::Map<const Eigen::VectorXd> ub_(ub.memptr(), ub.size());
  setAccelLimits(lb_, ub_);
}

void GMP_MPC::setAccelLimits(const Eigen::VectorXd &lb, const Eigen::VectorXd &ub)
{
  if (lb.size() != this->n_dof) throw std::runtime_error(GMP_MPC_fun_ + "lb.size() must be equal to n_dof");
  if (ub.size() != this->n_dof) throw std::runtime_error(GMP_MPC_fun_ + "ub.size() must be equal to n_dof");

  this->accel_lb = lb;
  this->accel_ub = ub;
}

// -----------------------------------------------

void GMP_MPC::setPosSlackLimit(double s_lim)
{   
  if (this->pos_slack) this->pos_slack_lim = s_lim*ones_ndof;
}

void GMP_MPC::setVelSlackLimit(double s_lim)
{    
  if (this->vel_slack) this->vel_slack_lim = s_lim*ones_ndof;
}

void GMP_MPC::setAccelSlackLimit(double s_lim)
{   
  if (this->accel_slack) this->accel_slack_lim = s_lim*ones_ndof;
}

void GMP_MPC::setCanonicalSystemFunction(std::function<std::array<double,2>(double, double)> can_sys_fun)
{
  this->can_sys_fun = can_sys_fun;
}

void GMP_MPC::integratePhase(double s, double s_dot, std::vector<double> &si_data, std::vector<double> &si_dot_data, std::vector<double> &si_ddot_data)
{
  si_data.resize(N);
  si_dot_data.resize(N);
  si_ddot_data.resize(N);

  std::vector<double> ti(this->N);
  ti[0] = 0;
  for (int i=1; i<ti.size(); i++) ti[i] = ti[i-1] + this->dt_[i-1];
  auto ode_fun = [this](const ode_state &state, ode_state &state_dot, double t)
                          { state_dot = this->can_sys_fun(state[0], state[1]); };
  controlled_stepper_type controlled_stepper;
  ode_state state = {s, s_dot};
  for (int i=0; i<ti.size()-1; i++)
  {
    si_data[i] = state[0];
    si_dot_data[i] = state[1];
    si_ddot_data[i] = ( this->can_sys_fun(state[0], state[1]) )[1];
    boost::numeric::odeint::integrate_adaptive(controlled_stepper, ode_fun, state, ti[i], ti[i+1], ti[i+1] - ti[i]);
  }
  si_data[this->N-1] = state[0];
  si_dot_data[this->N-1] = state[1];
  si_ddot_data[this->N-1] = ( this->can_sys_fun(state[0], state[1]) )[1];
}

GMP_MPC::Solution GMP_MPC::solve(double s, double s_dot)
{
  #ifdef GMP_MPC_L0G_STATISTICS
    timer.tic();
  #endif

  if (!this->can_sys_fun)
    throw std::runtime_error(GMP_MPC_fun_ + "The canonical system function has to be set at least once!");

  arma::vec y = arma::vec(&x0(0), n_dof, false);
  arma::vec y_dot = arma::vec(&x0(n_dof), n_dof, false);
  arma::vec y_ddot = arma::vec(&x0(2*n_dof), n_dof, false);
  arma::vec slack_var = arma::vec().zeros(n_slack);

  Eigen::VectorXd q = Eigen::VectorXd::Zero(n_var);

  Eigen::VectorXd lb_i(3*n_dof);
  lb_i << this->pos_lb, this->vel_lb, this->accel_lb;
  Eigen::VectorXd ub_i(3*n_dof);
  ub_i << this->pos_ub, this->vel_ub, this->accel_ub;

  Eigen::VectorXd slack_lim(n_slack);
  slack_lim << this->pos_slack_lim, this->vel_slack_lim, this->accel_slack_lim;

  // Eigen::SpMat Hi(n,n);
  std::vector<Eigen::Triplet<double>> H_triplets;
  H_triplets.reserve( N*n_var*n_var );

  Eigen::VectorXd qi = Eigen::VectorXd::Zero(n_var);

  std::vector<Eigen::Triplet<double>> Aineq_triplets;
  Aineq_triplets.reserve( n_ineq * n_var );

  //if (this->n_slack) Hi.submat(n-n_slack, n-n_slack, n-1, n-1) = this->Q_slack;

  // DMP phase variable
  std::vector<double> si_data, si_dot_data, si_ddot_data;
  integratePhase(s, s_dot, si_data, si_dot_data, si_ddot_data);
  double s_ddot = ( this->can_sys_fun(s, s_dot) )[1];

  #ifdef GMP_MPC_L0G_STATISTICS
    timer_H.tic();
  #endif

  std::vector<Eigen::Triplet<double>> A_obst_triplets;
  A_obst_triplets.reserve( n_obst * n_var );
  Eigen::VectorXd b_obst(n_obst);
  Eigen::VectorXd b_obst_up = 1e30*Eigen::VectorXd::Ones(n_obst);

  log.clear();

  for (int i=0; i<this->N; i++)
  {
    double si = si_data[i];
    double si_dot = si_dot_data[i];
    double si_ddot = si_ddot_data[i];

    arma::vec yd_i_ = this->gmp_ref->getYd(si);
    arma::vec dyd_i_ = this->gmp_ref->getYdDot(si, si_dot);
    Eigen::Map<Eigen::VectorXd> yd_i(yd_i_.memptr(), n_dof);
    Eigen::Map<Eigen::VectorXd> dyd_i(dyd_i_.memptr(), n_dof);

    if (si >= 1)
    {
      yd_i = x_f.segment(0,n_dof);
      dyd_i = x_f.segment(n_dof, n_dof);
    }

    arma::vec phi = this->gmp_mpc->regressVec(si);
    arma::vec phi_dot = this->gmp_mpc->regressVecDot(si, si_dot);
    arma::vec phi_ddot = this->gmp_mpc->regressVecDDot(si, si_dot, si_ddot);

    // Eigen::SpMat Qi_;
    // if (i==this->N-1) Qi_ = this->QN;
    // else Qi_ = this->Qi;

    Eigen::SpMat Psi_gain = createPsiMat(n_dof, sqrt_pos_gain*phi, sqrt_vel_gain*phi_dot);
    Eigen::VectorXd xd_i(2*n_dof);
    xd_i << yd_i, dyd_i;

    // Eigen::SpMat Hi_ = Psi.transpose()*Qi_*Psi;
    Eigen::SpMat Hi_ = Psi_gain.transpose()*Psi_gain;

    Eigen::appendTriplets(H_triplets, Hi_);
    // Hi.submat(0,0, n-n_slack-1, n-n_slack-1) =
    // H = H + Hi;

    // std::cout << "Hi_:\n" << std::setprecision(2) << std::setw(4) << Hi_.toDense() << "\n";

    //qi.segment(0, n_var-n_slack) = -Psi.transpose()*Qi_*xd_i;
    qi.segment(0, n_var-n_slack) = -Psi_gain.transpose()*(xd_i.cwiseProduct(sqrt_pos_vel_gain_vec));
    q = q + qi;

    // arma::mat Aineq_i = arma::join_vert( arma::kron(I_ndof,phi.t()), arma::kron(I_ndof,phi_dot.t()), arma::kron(I_ndof,phi_ddot.t()) );
    // Aineq.rows(i*n_dof3, (i+1)*n_dof3-1) = arma::join_horiz(Aineq_i, this->Aineq_slack);
    addToAineqTriplets(Aineq_triplets, n_dof, phi, phi_dot, phi_ddot, i*n_dof3);

    // triplets for identity to insert slacks in constraints: [Aj I] * [w; s]
    addTriplets(Aineq_triplets, Aineq_slack_triplets, i*n_dof3, n_var-n_slack);

    // check collisions with ellipsoids and calculate plane constraints to avoid them
    arma::rowvec Ai;
    double bi;
    arma::vec y_pred = this->getYd(si);
    this->process_obstacles(y_pred, phi, &Ai, &bi);
    addTriplets(A_obst_triplets, Ai, i, 0);
    b_obst(i) = bi;

    log.yd_points.emplace_back(arma::vec(yd_i.data(), yd_i.size(), true));
    log.y_pred_points.push_back(y_pred);

    // si = si + si_dot*this->dt_(i-1);
    // si_dot = si_dot + si_ddot*this->dt_(i-1);
    // si_ddot = ... (if it changes too)
  }

  log.y_current = {x0(0), x0(1), x0(2)};
  log.yd_current = this->gmp_ref->getYd(s);
  log.si_data = si_data;
  log.W_opt = W_mpc;
  log.Wd = this->gmp_ref->W;
  log.y_target = arma::vec(x_f.segment(0,n_dof).data(), n_dof, false);
  if (plot_callback) plot_callback(log);

  // =======  initial state constraint  =======
  Eigen::SpMat Aeq = Eigen::join_horiz( this->Phi0, Eigen::SpMat(n_dof3, this->n_slack) );
  Eigen::VectorXd beq = this->x0;

  // =======  final state constraint  =======
  // Aineq = arma::join_vert(Aineq, arma::join_horiz(this->Phi_f, arma::mat().zeros(n_dof3, this->n_slack) ) );
  addTriplets(Aineq_triplets, Phi_f_triplets, N*n_dof3, 0);
  // addToAineqTriplets(Aineq_triplets, n_dof, phi_f, phi_f_dot, phi_f_ddot, N*n_dof3);
  Eigen::VectorXd lb_f = this->x_f - this->err_tol_f;
  Eigen::VectorXd ub_f = this->x_f + this->err_tol_f;

  // add final state to cost function too, to help the optimizer
  // Maybe no... it could affect tracking during intermediate steps...
//  Eigen::SpMat Hi_ = Phi_f.transpose()*Qf*Phi_f;
//  Eigen::appendTriplets(H_triplets, Hi_);
//  q.segment(0, n_var - n_slack) -= Phi_f.transpose()*Qf*this->x_f;

  // =======  process via-points  =======

  // discard past (non active) via-points
  for (auto it = via_points.begin(); it != via_points.end(); )
  {
    if ( it->s < s ) it = via_points.erase(it);
    else it++;
  }

  // consider only up to the first 'n_via' active via-points
  unsigned n_vp = std::min( (size_t)n_via, this->via_points.size() );

  Eigen::VectorXd lb_v = -inf*Eigen::VectorXd::Ones(n_via * n_dof3);
  Eigen::VectorXd  ub_v = -lb_v;

  auto vp_it = via_points.begin();
  for (int i=0; i<n_vp; i++)
  {
      ViaPoint vp = *vp_it++;
      arma::vec phi = this->gmp_mpc->regressVec(vp.s);
      arma::vec phi_dot = this->gmp_mpc->regressVecDot(vp.s, s_dot);
      arma::vec phi_ddot = this->gmp_mpc->regressVecDDot(vp.s, s_dot, 0); // ? s_ddot

      addToAineqTriplets(Aineq_triplets, n_dof, phi, phi_dot, phi_ddot, (N+1+i)*n_dof3);
      lb_v.segment(i*n_dof3, n_dof3) << (vp.pos.array() - vp.err_tol), vel_lb, accel_lb;
      ub_v.segment(i*n_dof3, n_dof3) << (vp.pos.array() + vp.err_tol), vel_ub, accel_ub;

      // add via-point to cost function too, to help the optimizer
      double Qi_ = 1e3;
      Eigen::appendTriplets(H_triplets, vp.Phi.transpose()*(Qi_*vp.Phi) );
      q.segment(0, n_var - n_slack) -= vp.Phi.transpose()*(Qi_*vp.pos);
  }

  // // =======  expand/shrink dual solution guess, in case constraints were added/removed =======
  // int n_ineq_plus = Aineq.n_rows - this->Z0_dual_ineq.size();
  // if (n_ineq_plus>0) this->Z0_dual_ineq = arma::join_vert(this->Z0_dual_ineq, arma::vec().zeros(n_ineq_plus) );
  // else this->Z0_dual_ineq.resize(Aineq.n_rows);

  // add diagonal to H for numerical stability
  addDiagonalTriplets(H_triplets, n_var, 1e-6);
  addDiagonalTriplets(H_triplets, N*Qs_vec, n_var-n_slack, n_var-n_slack);
  Eigen::SpMat H(n_var, n_var);
  H.setFromTriplets(H_triplets.begin(), H_triplets.end());
  // H.setIdentity();
  // H.diagonal() *= 1e-6; // for numerical stability

  // =======  obstacle constraints triplets  =======
  addTriplets(Aineq_triplets, A_obst_triplets, n_ineq-n_obst-n_slack, 0);

//  Eigen::SpMat A_obst(n_obst, n_var);
//  A_obst.setFromTriplets(A_obst_triplets.begin(), A_obst_triplets.end());
//  std::cout << "A_obst:\n" << A_obst.toDense() << "\nb_obst: " << b_obst.transpose() << "\n";

  // triplets for slack variable limits
  addDiagonalTriplets(Aineq_triplets, n_slack, 1.0, n_ineq-n_slack, n_var-n_slack);
  // if (this->n_slack)
  // {
  //   int i_end = Aineq.rows() - 1;
  //   int j_end = Aineq.cols() - 1;
  //   // identity matrix for slack variable bounds
  //   Aineq.submat(i_end-this->n_slack+1, j_end-this->n_slack+1, i_end, j_end) = arma::mat().eye(this->n_slack,this->n_slack);
  // }

  Eigen::SpMat Aineq(n_ineq, n_var);
  Aineq.setFromTriplets(Aineq_triplets.begin(), Aineq_triplets.end());

  Eigen::VectorXd lb(n_ineq);
  lb << lb_i.replicate(this->N,1), lb_f, lb_v, b_obst, -slack_lim;
  Eigen::VectorXd ub(n_ineq);
  ub << ub_i.replicate(this->N,1), ub_f, ub_v, b_obst_up, slack_lim;

//  Eigen::MatrixXd lb_A_ub( lb.rows() , lb.cols() + Aineq.cols() + ub.cols() );
//  lb_A_ub << lb, Aineq.toDense(), ub;
//  Eigen::MatrixXd Aeq_beq( beq.rows() , beq.cols() + Aeq.cols() );
//  Aeq_beq << Aeq.toDense(), beq;
//
//   std::cout << "========== QP problem ============\n";
//   std::cout << std::setprecision(2) << std::setw(4)
//             << "H: \n" << H.toDense() << "\n"
//             << "q: \n" << q.transpose() << "\n"
//             << "lb < A < ub :\n" << lb_A_ub << "\n"
//             << "Aeq = beq: \n" << Aeq_beq << "\n";

  // exit(-1);

  #ifdef GMP_MPC_L0G_STATISTICS
    calcH_elaps_t_data.push_back(timer_H.toc() * 1000);
  #endif

  #ifdef GMP_MPC_L0G_STATISTICS
    if (s > s_rec)
    {
      s_rec = 1e10;
      this->H_ = H;
      this->A_ = Aineq;
    }
  #endif

  // ===========  solve optimization problem  ==========
  #ifdef GMP_MPC_L0G_STATISTICS
    timer_solve.tic();
  #endif

  osqp_::QPsolver solver(H, q, Aineq, lb, ub, Aeq, beq);
  solver.setPrimalSolutionGuess(this->Z0);
  solver.setDualSolutionGuess(this->Z0_dual_ineq, this->Z0_dual_eq);

  solver.settings->time_limit = this->settings.time_limit;
  solver.settings->max_iter = this->settings.max_iter;
  solver.settings->eps_abs = this->settings.abs_tol;
  solver.settings->eps_rel = this->settings.rel_tol;
  solver.settings->warm_start = false;
  solver.settings->verbose = false;

  int exit_flag = solver.solve();

  #ifdef GMP_MPC_L0G_STATISTICS
    solve_elaps_t_data.push_back( timer_solve.toc()*1000 );
  #endif

  std::string exit_msg = solver.getExitMsg();
  if (exit_flag >= 0)
  {
    // =======  extract the solution  =======
    arma::vec X = solver.getPrimalSolution();
    arma::vec w = X.subvec(0, X.n_elem-1-this->n_slack);
    if (this->n_slack) slack_var = X.subvec(X.n_elem-this->n_slack, X.n_elem-1);
    this->W_mpc = arma::reshape(w, N_kernels, n_dof).t();

    this->Z0 = X;
    this->Z0_dual_ineq = solver.getIneqDualSolution();
    this->Z0_dual_eq = solver.getEqDualSolution();

    // =======  Generate optimal output  =======
    y = this->getYd(s);
    y_dot = this->getYdDot(s, s_dot);
    y_ddot = this->getYdDDot(s, s_dot, s_ddot);

    if (this->n_slack)
    {
      // check, because if the optimization finishes prematurely,
      // sometimes it can happen that the slacks are violated...
      // Not sure why this happens...?
      double max_violation = arma::max( arma::abs(slack_var) - arma::vec(slack_lim.data(), slack_lim.size(), false) );
      if ( max_violation > 5e-4 )
      {
        exit_msg = "Slack variable limits violated: max violation = " + std::to_string(max_violation);
        exit_flag = -1;
        goto end_solve;
      }
    }

    this->setInitialState(y, y_dot, y_ddot, s, s_dot, s_ddot);
  }

  #ifdef GMP_MPC_L0G_STATISTICS
    elaps_t_data.push_back( timer.toc()*1000 );
  #endif

end_solve:
  Solution sol;
  sol.y = y;
  sol.y_dot = y_dot;
  sol.y_ddot = y_ddot;
  int i1=0, i2=n_dof-1;
  if (pos_slack)
  {
    sol.pos_slack = slack_var.subvec(i1,i2);
    i1 += n_dof;
    i2 += n_dof;
  }
  if (vel_slack)
  {
    sol.vel_slack = slack_var.subvec(i1,i2);
    i1 += n_dof;
    i2 += n_dof;
  }
  if (accel_slack) sol.accel_slack = slack_var.subvec(i1,i2);
  sol.exit_flag = exit_flag;
  sol.exit_msg = exit_msg;
  return sol;
}

arma::vec GMP_MPC::getYd(double s) const
{
  return W_mpc*gmp_mpc->regressVec(s);
}

arma::vec GMP_MPC::getYdDot(double s, double s_dot) const
{
  return W_mpc*gmp_mpc->regressVecDot(s, s_dot);
}

arma::vec GMP_MPC::getYdDDot(double s, double s_dot, double s_ddot) const
{
  return W_mpc*gmp_mpc->regressVecDDot(s, s_dot, s_ddot);
}
  
void GMP_MPC::process_obstacles(const arma::vec &p1, const arma::vec &phi, arma::rowvec *Ai, double *bi)
{
  // check collisions with ellipsoids and calculate plane constraints to avoid them
  Ai->zeros(this->n_var);
  *bi = -1;
  for (auto &it : this->obst_map)
  {
    auto obst = it.second;
    if (obst->getConstrPlain(p1, phi, Ai, bi))
    {
      *Ai = arma::join_horiz(*Ai, arma::rowvec().zeros(this->n_slack));

      log.n_e_data.push_back(obst->n_e);
      log.p_e_data.push_back(obst->p_e);
      log.p_data.push_back(obst->p);

      break;
    }
  }
}

} // namespace gmp_

} // namespace as64_
