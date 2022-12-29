#ifndef OSQP_QUADPROG_H
# define OSQP_QUADPROG_H

#include <osqp_lib/osqp.h>
#include <osqp_lib/csc_mat.h>
#include <armadillo>

#ifndef NULL_MAT
  #define NULL_MAT arma::mat(0,0)
#endif

namespace osqp_
{

enum QuadProgSolutionStatus
{
  OPTIMAL = 0,
  SUBOPTIMAL,
  FAILED
};

class QuadProgOptions 
{
public:
  QuadProgOptions()
  {
    max_iters = 2000;
    warm_start = true;
    polish = false;
    time_limit = 0;
    verbose = false;
    parallel = false;
  }

  int max_iters;
  bool warm_start;
  bool polish;
  double time_limit;
  bool verbose;
  bool parallel;
};

typedef struct
{
  arma::mat x;
  QuadProgSolutionStatus status;
  std::string exit_msg;

} QuadProgSolution;

QuadProgSolution quadprog(const arma::mat &H, const arma::mat &f, const arma::mat &A=NULL_MAT, const arma::mat &lb=NULL_MAT, const arma::mat &ub=NULL_MAT,
              const arma::mat &Aeq=NULL_MAT, const arma::mat &beq=NULL_MAT, const QuadProgOptions &options=QuadProgOptions());

} // namespace osqp_

#endif // OSQP_QUADPROG_H
