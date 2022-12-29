
#include <osqp_lib/csc_mat.h>

namespace osqp_
{

CSC_mat::CSC_mat(const arma::mat &A, bool upper_diag_mat)
{
  if (upper_diag_mat) process_UD_mat(A);
  else process_mat(A);
}

void CSC_mat::process_mat(const arma::mat &A)
{
  int m = A.n_rows;
  int n = A.n_cols;

  n_rows = m;
  n_cols = n;

  data.resize(m*n);
  row_ind.resize(m*n);
  cs.resize(n+1);

  cs[0] = 0;

  nnz = 0;

  int k=0;

  for (int j=0; j<n; j++)
  {
    int s = 0;
    for (int i=0; i<m ;i++)
    {
      if (A(i,j))
      {
        data[k] = A(i,j);
        row_ind[k] = i;
        k++;
        s++;
      }
    }
    nnz += s;
    cs[j+1] = nnz; // cs[j] + s;
  }

  data.resize(k);
  row_ind.resize(k);
}

void CSC_mat::process_UD_mat(const arma::mat &A)
{
  int m = A.n_rows;
  int n = A.n_cols;

  if (m != n) throw std::runtime_error("[CSC_mat::process_UD_mat]: The matrix must be square.");

  n_rows = m;
  n_cols = n;

  int n_max = (n+1)*n/2;

  data.resize(n_max);
  row_ind.resize(n_max);
  cs.resize(n+1);

  cs[0] = 0;

  nnz = 0;

  int k=0;

  for (int j=0; j<n; j++)
  {
    int s = 0;
    for (int i=0; i<=j ;i++)
    {
      if (A(i,j))
      {
        data[k] = A(i,j);
        row_ind[k] = i;
        k++;
        s++;
      }
    }
    nnz += s;
    cs[j+1] = nnz; // cs[j] + s;
  }

  data.resize(k);
  row_ind.resize(k);
}

} // namespace osqp_

