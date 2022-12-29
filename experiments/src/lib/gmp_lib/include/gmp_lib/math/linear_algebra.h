#ifndef GMP_LIB_LIN_ALG__64_H
#define GMP_LIB_LIN_ALG__64_H

#include <armadillo>

#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace Eigen
{

typedef SparseMatrix<double, ColMajor> SpMat;
typedef SparseVector<double> SpVec;

inline void appendTriplets(std::vector<Eigen::Triplet<double>> &values, const SpMat &sm, int row_offset=0, int col_offset=0)
{
  for (int k=0; k<sm.outerSize(); ++k)
  {
    for (SparseMatrix<double>::InnerIterator it(sm,k); it; ++it)
      values.push_back( Eigen::Triplet<double>( it.row()+row_offset, it.col()+col_offset, it.value()) );
  }
}


inline SpMat blkdiag(const SpMat &s1, const SpMat &s2)
{
  int m1 = s1.rows();
  int n1 = s1.cols();

  std::vector<Eigen::Triplet<double>> values;
  values.reserve( s1.nonZeros() + s2.nonZeros() );
  appendTriplets(values, s1);
  appendTriplets(values, s2, m1, n1);
  
  SpMat s(m1+s2.rows(), n1+s2.cols());
  s.setFromTriplets(values.begin(), values.end());
  return s;
}

inline SpMat blkdiag(const SpMat &s1, const SpMat &s2, const SpMat &s3)
{
  int m1 = s1.rows();
  int n1 = s1.cols();

  int m2 = s2.rows();
  int n2 = s2.cols();

  std::vector<Eigen::Triplet<double>> values;
  values.reserve( s1.nonZeros() + s2.nonZeros() + s3.nonZeros() );
  appendTriplets(values, s1);
  appendTriplets(values, s2, m1, n1);
  appendTriplets(values, s3, m1+m2, n1+n2);
  
  SpMat s(m1+m2+s3.rows(), n1+n2+s3.cols());
  s.setFromTriplets(values.begin(), values.end());
  return s;
}

inline SpMat join_horiz(const SpMat &s1, const SpMat &s2)
{
  if (s1.rows() != s2.rows())
    throw std::runtime_error("[Eigen::join_horiz]: Dimensions of matrices are not consistent.");

  int m = s1.rows();
  int n1 = s1.cols();

  std::vector<Eigen::Triplet<double>> values;
  values.reserve( s1.nonZeros() + s2.nonZeros() );
  appendTriplets(values, s1);
  appendTriplets(values, s2, 0, n1);
  
  SpMat s(m, n1+s2.cols());
  s.setFromTriplets(values.begin(), values.end());
  return s;
}

inline SpMat join_horiz(const SpMat &s1, const SpMat &s2, const SpMat &s3)
{
  if (s1.rows() != s2.rows() || s1.rows() != s3.rows())
    throw std::runtime_error("[Eigen::join_horiz]: Dimensions of matrices are not consistent.");

  int m = s1.rows();

  int n1 = s1.cols();
  int n2 = s2.cols();

  std::vector<Eigen::Triplet<double>> values;
  values.reserve( s1.nonZeros() + s2.nonZeros() + s3.nonZeros() );
  appendTriplets(values, s1);
  appendTriplets(values, s2, 0, n1);
  appendTriplets(values, s3, 0, n1+n2);
  
  SpMat s(m, n1+n2+s3.cols());
  s.setFromTriplets(values.begin(), values.end());
  return s;
}

inline SpMat join_vert(const SpMat &s1, const SpMat &s2)
{ 
  if (s1.cols() != s2.cols())
    throw std::runtime_error("[Eigen::join_vert]: Dimensions of matrices are not consistent.");

  int n = s1.cols();
  int m1 = s1.rows();

  std::vector<Eigen::Triplet<double>> values;
  values.reserve( s1.nonZeros() + s2.nonZeros() );
  appendTriplets(values, s1);
  appendTriplets(values, s2, m1, 0);
  
  SpMat s(m1+s2.rows(), n);
  s.setFromTriplets(values.begin(), values.end());
  return s;
}

inline SpMat join_vert(const SpMat &s1, const SpMat &s2, const SpMat &s3)
{
  if (s1.cols() != s2.cols() || s1.cols() != s3.cols())
    throw std::runtime_error("[Eigen::join_vert]: Dimensions of matrices are not consistent.");

  int n = s1.cols();

  int m1 = s1.rows();
  int m2 = s2.rows();

  std::vector<Eigen::Triplet<double>> values;
  values.reserve( s1.nonZeros() + s2.nonZeros() + s3.nonZeros() );
  appendTriplets(values, s1);
  appendTriplets(values, s2, m1, 0);
  appendTriplets(values, s3, m1+m2, 0);
  
  SpMat s(m1+m2+s3.rows(), n);
  s.setFromTriplets(values.begin(), values.end());
  return s;
}

// template<typename T>
// Matrix<T, Dynamic, Dynamic> join_horiz(const Matrix<T, Dynamic, Dynamic> &s1, const Matrix<T, Dynamic, Dynamic> &s2)
// {
//   if (s1.rows() != s2.rows())
//     throw std::runtime_error("[Eigen::join_horiz]: Dimensions of matrices are not consistent.");

//   Matrix<T, Dynamic, Dynamic> s( s1.rows(), s1.cols() + s2.cols() );
//   s << s1, s2;
//   return s;
// }


// template<typename T>
// Matrix<T, Dynamic, Dynamic> join_vert(const Matrix<T, Dynamic, Dynamic> &s1, const Matrix<T, Dynamic, Dynamic> &s2)
// {
//   if (s1.cols() != s2.cols())
//     throw std::runtime_error("[Eigen::join_vert]: Dimensions of matrices are not consistent.");

//   Matrix<T, Dynamic, Dynamic> s( s1.rows() + s2.rows(), s1.cols() );
//   s << s1, s2;
//   return s;
// }

// template<typename T>
// Matrix<T, Dynamic, Dynamic> join_vert(const Matrix<T, Dynamic, Dynamic> &s1, const Matrix<T, Dynamic, Dynamic> &s2)
// {


// }

// template<typename T>
// Matrix<T, Dynamic, 1> join_vert(const Matrix<T, Dynamic, 1> &s1, const Matrix<T, Dynamic, 1> &s2)
// {
//   Matrix<T, Dynamic, 1> s( s1.rows() + s2.rows() );
//   s << s1, s2;
//   return s;
// }

template<typename Derived>
EigenBase<Derived> join_horiz(const EigenBase<Derived>& m1, const EigenBase<Derived>& m2)
{
  if (m1.rows() != m2.cols()) throw std::runtime_error("[Eigen::join_horiz]:Rows of matrices being concatenated are inconsistent: " + std::to_string(m1.rows()) + std::to_string(m2.rows()) );
  
  EigenBase<Derived> m( m1.rows()+m2.rows() , m1.cols()+m2.cols() );
  m << m1, m2;
  return m;
}

} // namespace Eigen


#endif // GMP_LIB_LIN_ALG__64_H