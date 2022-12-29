#ifndef OSQP_QP_SOLVER_H
# define OSQP_QP_SOLVER_H

#include <memory>
#include <osqp_lib/osqp.h>

#include <armadillo>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Eigen
{
  typedef Eigen::SparseMatrix<double, ColMajor> SpMat;
} // namespace Eigen

#ifndef NULL_MAT
  #define NULL_MAT arma::mat(0,0)
#endif

#ifndef NULL_VEC
  #define NULL_VEC arma::vec(0)
#endif

namespace osqp_
{

class QPsolver
{
public:
  QPsolver(const arma::mat &H, const arma::vec &q,
           const arma::mat &A, const arma::vec &lb, const arma::vec &ub,
           const arma::mat &Aeq, const arma::vec &beq);

  QPsolver(const Eigen::SpMat &H, const Eigen::VectorXd &q,
           const Eigen::SpMat &A, const Eigen::VectorXd &lb, const Eigen::VectorXd &ub,
           const Eigen::SpMat &Aeq, const Eigen::VectorXd &beq);

  ~QPsolver();

  void setPrimalSolutionGuess(const arma::vec &x0);
  void setPrimalSolutionGuess(const Eigen::VectorXd &x0);

  void setDualSolutionGuess(const arma::vec &y0_ineq, const arma::vec &y0_eq);
  void setDualSolutionGuess(const Eigen::VectorXd &y0_ineq, const Eigen::VectorXd &y0_eq);

  int solve();

  arma::vec getPrimalSolution() const;

  arma::vec getIneqDualSolution() const;

  arma::vec getEqDualSolution() const;

  std::string getExitMsg() const;

  OSQPSettings  *settings;

private:

  void init(const Eigen::SpMat &H, const Eigen::VectorXd &q,
           const Eigen::SpMat &A, const Eigen::VectorXd &lb, const Eigen::VectorXd &ub,
           const Eigen::SpMat &Aeq, const Eigen::VectorXd &beq);

  struct CSC_struct
  {
    void fill(const Eigen::SpMat &sm)
    {
      int nnz = sm.nonZeros();
      int n_cols = sm.cols();

      values.reserve(nnz);
      inner_ind.reserve(nnz);
      outer_ind.resize(n_cols+1);

      const double *value_ptr = sm.valuePtr();
      const int *inner_ind_ptr = sm.innerIndexPtr();
      const int *outer_ind_ptr = sm.outerIndexPtr();
      const int *inner_nnz_ptr = sm.innerNonZeroPtr();

      for (int j=0; j<n_cols; j++)
      {
        outer_ind[j] = outer_ind_ptr[j];
        int k1 = outer_ind_ptr[j];
        int k2 = outer_ind_ptr[j+1];
        if (inner_nnz_ptr) k2 = k1 + inner_nnz_ptr[j];

        for (int k=k1; k<k2; k++)
        {
          values.push_back( value_ptr[k] );
          inner_ind.push_back( inner_ind_ptr[k] );
        }
      }
      outer_ind[n_cols] = outer_ind_ptr[n_cols];;

    }

    unsigned nnz() { return values.size(); }
    c_float *valuesPtr() { return values.data(); }
    c_int *innerIndPtr() { return inner_ind.data(); }
    c_int *outerIndPtr() { return outer_ind.data(); }

    std::vector<c_float> values;
    std::vector<c_int> inner_ind;
    std::vector<c_int> outer_ind;
  };

  CSC_struct P_csc;
  CSC_struct A_csc;

  Eigen::VectorXd lb_;
  Eigen::VectorXd ub_;
  Eigen::VectorXd q_;

  OSQPWorkspace *work;
  OSQPData      *data;

  int n_ineq;
  int n_eq;

  bool warm_start_x;
  arma::vec x0;

  bool warm_start_y;
  arma::vec y0_ineq, y0_eq;

  std::string exit_msg;

}; // class QPsolver


} // namespace osqp_

#endif // OSQP_QP_SOLVER_H
