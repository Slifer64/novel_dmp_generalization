#ifndef OSQP_CSC_MAT_H
# define OSQP_CSC_MAT_H

#include <osqp_lib/glob_opts.h>
#include <armadillo>

namespace osqp_
{

class CSC_mat
{
public:
  CSC_mat() {}
  CSC_mat(const arma::mat &A, bool upper_diag_mat=false);

  int n_rows;  // number of rows
  int n_cols;  // number of comumns
  int nnz; // number of nonzeros
  std::vector<c_float> data;  // value of each non-zero element
  std::vector<c_int> row_ind;  // row index of each non-zero element
  std::vector<c_int> cs;  // cumulative sum of none-zero elements from column to column

  c_float *dataPtr() { return &data[0]; }
  c_int *rowIndPtr() { return &row_ind[0]; }
  c_int *csPtr() { return &cs[0]; }

private:
  void process_mat(const arma::mat &A);
  void process_UD_mat(const arma::mat &A);
}; // CSC_mat

} // namespace osqp_

#endif // OSQP_CSC_MAT_H
