#include <osqp_lib/qp_solver.h>
#include <osqp_lib/csc_mat.h>


namespace osqp_
{

#define QPsolver_fun_ std::string("[QPsolver::") + __func__ + "]: "

Eigen::SpMat join_vert(const Eigen::SpMat &s1, const Eigen::SpMat &s2)
{ 
  if (s1.cols() != s2.cols())
    throw std::runtime_error("[OSQP::join_vert]: Dimensions of matrices are not consistent.");

  int n = s1.cols();
  int m1 = s1.rows();

  std::vector<Eigen::Triplet<double>> values;
  values.reserve( s1.nonZeros() + s2.nonZeros() );

  for (int k=0; k<s1.outerSize(); ++k)
  {
    for (Eigen::SpMat::InnerIterator it(s1,k); it; ++it)
      values.push_back( Eigen::Triplet<double>( it.row(), it.col(), it.value()) );
  }
  for (int k=0; k<s2.outerSize(); ++k)
  {
    for (Eigen::SpMat::InnerIterator it(s2,k); it; ++it)
      values.push_back( Eigen::Triplet<double>( it.row()+m1, it.col(), it.value()) );
  }
  
  Eigen::SpMat s(m1+s2.rows(), n);
  s.setFromTriplets(values.begin(), values.end());
  return s;
}

QPsolver::QPsolver(const arma::mat &H, const arma::vec &q,
         const arma::mat &A, const arma::vec &lb, const arma::vec &ub,
         const arma::mat &Aeq, const arma::vec &beq)
{

  Eigen::SpMat H_ = Eigen::Map<const Eigen::MatrixXd>(H.memptr(), H.n_rows, H.n_cols).sparseView();
  Eigen::SpMat A_ = Eigen::Map<const Eigen::MatrixXd>(A.memptr(), A.n_rows, A.n_cols).sparseView();
  Eigen::SpMat Aeq_ = Eigen::Map<const Eigen::MatrixXd>(Aeq.memptr(), Aeq.n_rows, Aeq.n_cols).sparseView();

  Eigen::Map<const Eigen::VectorXd> q_(q.memptr(), q.size());
  Eigen::Map<const Eigen::VectorXd> lb_(lb.memptr(), lb.size());
  Eigen::Map<const Eigen::VectorXd> ub_(ub.memptr(), ub.size());
  Eigen::Map<const Eigen::VectorXd> beq_(beq.memptr(), beq.size());

  init(H_, q_, A_, lb_, ub_, Aeq_, beq_);
}

QPsolver::QPsolver(const Eigen::SpMat &H, const Eigen::VectorXd &q,
           const Eigen::SpMat &A, const Eigen::VectorXd &lb, const Eigen::VectorXd &ub,
           const Eigen::SpMat &Aeq, const Eigen::VectorXd &beq)
{
  init(H, q, A, lb, ub, Aeq, beq);
}

void QPsolver::init(const Eigen::SpMat &H, const Eigen::VectorXd &q,
           const Eigen::SpMat &A, const Eigen::VectorXd &lb, const Eigen::VectorXd &ub,
           const Eigen::SpMat &Aeq, const Eigen::VectorXd &beq)
{
  work = NULL;
  settings = NULL;
  data = NULL;
  warm_start_x = false;
  warm_start_y = false;

  if (H.rows() != H.cols()) throw std::runtime_error(QPsolver_fun_ + "H must be square");
  if (H.rows() != q.size()) throw std::runtime_error(QPsolver_fun_ + "H.n_cols must be equal to q.size()");
  if (A.size())
  {
    if (A.cols() != q.size()) throw std::runtime_error(QPsolver_fun_ + "A.n_cols must be equal to q.size()");
    if (A.rows() != lb.size()) throw std::runtime_error(QPsolver_fun_ + "A.n_rows must be equal to lb.size()");
    if (A.rows() != ub.size()) throw std::runtime_error(QPsolver_fun_ + "A.n_rows must be equal to ub.size()");
  }
  if (Aeq.size())
  {
    if (Aeq.cols() != q.size()) throw std::runtime_error(QPsolver_fun_ + "Aeq.n_cols must be equal to q.size()");
    if (Aeq.rows() != beq.size()) throw std::runtime_error(QPsolver_fun_ + "Aeq.n_rows must be equal to beq.size()");
  }

  Eigen::SpMat H_ = H.triangularView<Eigen::Upper>();
  P_csc.fill(H_);

  this->n_ineq = A.rows();
  this->n_eq = Aeq.rows();

  int n_var = q.size();
  int n_constr = A.rows() + Aeq.rows();

  Eigen::SpMat Aineq = join_vert(A, Aeq);
  A_csc.fill(Aineq);

  lb_.resize(n_constr);
  lb_ << lb, beq;

  ub_.resize(n_constr);
  ub_ << ub, beq;

  q_ = q;

  settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
  if (!settings) throw std::runtime_error(QPsolver_fun_ + "Failed to allocate space for OSQPSettings...");

  data = (OSQPData *)c_malloc(sizeof(OSQPData));
  if (!data) throw std::runtime_error(QPsolver_fun_ + "Failed to allocate space for QSQPData...");

  data->n = n_var;
  data->m = n_constr;
  data->P = csc_matrix(data->n, data->n, P_csc.nnz(), P_csc.valuesPtr(), P_csc.innerIndPtr(), P_csc.outerIndPtr() );
  data->A = csc_matrix(data->m, data->n, A_csc.nnz(), A_csc.valuesPtr(), A_csc.innerIndPtr(), A_csc.outerIndPtr() );
  data->q = q_.data();
  data->l = lb_.data();
  data->u = ub_.data();

  osqp_set_default_settings(settings);
  // settings->alpha = 1.0; // Change alpha parameter
  settings->warm_start = false;
  settings->polish = 0;
  //settings->time_limit = 0;
  //settings->max_iter = 4000;
  settings->verbose = false;
  settings->warm_start = false;
}

QPsolver::~QPsolver()
{
  if (data)
  {
    if (data->A) c_free(data->A);
    if (data->P) c_free(data->P);
    c_free(data);
  }
  if (settings) c_free(settings);

  if (work) osqp_cleanup(work);
}

void QPsolver::setPrimalSolutionGuess(const arma::vec &x0)
{
  warm_start_x = true;
  this->x0 = x0;
}

void QPsolver::setPrimalSolutionGuess(const Eigen::VectorXd &x)
{
  setPrimalSolutionGuess( arma::vec(const_cast<double *>(x.data()), x.size(), false) );
}

void QPsolver::setDualSolutionGuess(const arma::vec &y0_ineq, const arma::vec &y0_eq)
{
  warm_start_y = true;
  this->y0_ineq = y0_ineq;
  this->y0_eq = y0_eq;
}

void QPsolver::setDualSolutionGuess(const Eigen::VectorXd &y0_ineq, const Eigen::VectorXd &y0_eq)
{
  setDualSolutionGuess( arma::vec(const_cast<double *>(y0_ineq.data()), y0_ineq.size(), false) , 
                        arma::vec(const_cast<double *>(y0_eq.data()), y0_eq.size(), false) );
}

int QPsolver::solve()
{
  c_int exit_flag = osqp_setup(&work, data, settings);
  if (exit_flag != 0)
  {
    exit_msg = work->info->status;
    return -1;
  }

  arma::vec y0 = arma::join_vert(y0_ineq, y0_eq);
  if (warm_start_x && warm_start_y) exit_flag = osqp_warm_start(work, x0.memptr(), y0.memptr());
  else if (warm_start_x) exit_flag = osqp_warm_start_x(work, x0.memptr());
  else if (warm_start_y) exit_flag = osqp_warm_start_y(work, y0.memptr());


  if (exit_flag != 0)
  {
    exit_msg = work->info->status;
    return -1;
  }

  // Solve Problem
  osqp_solve(work);

  exit_msg = work->info->status;

  c_int sol_status = work->info->status_val;

  if (sol_status == -3 || sol_status == -4 || sol_status == -7 || sol_status == -10) return -1;

  x0 = arma::mat(work->solution->x, data->n, 1, true);
  y0_ineq = arma::mat(work->solution->y, n_ineq, 1, true);
  y0_eq = arma::mat(work->solution->y+n_ineq, n_eq, 1, true);

  if (sol_status == 1) return 0;
  else return 1;
}

arma::vec QPsolver::getPrimalSolution() const
{
  return x0;
}

arma::vec QPsolver::getIneqDualSolution() const
{
  return y0_ineq;
}

arma::vec QPsolver::getEqDualSolution() const
{
  return y0_eq;
}

std::string QPsolver::getExitMsg() const
{
  return exit_msg;
}

} // namespace osqp_

