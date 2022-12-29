#include <gmp_lib/math/math.h>

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

namespace as64_
{

namespace gmp_
{

// =======================================================
// =======================================================

Eigen::Vector4d rotm2quat(Eigen::Matrix3d rotm, bool is_rotm_orthonormal)
{
    if (!is_rotm_orthonormal) return mat2quat(rotm);

    Eigen::Quaternion<double> temp_quat(rotm);
    Eigen::Vector4d quat;
    quat << temp_quat.w(), temp_quat.x(), temp_quat.y(), temp_quat.z();

    quat = quat * (2*(quat(0)>=0)-1); // to avoid discontinuities

    return quat;
}

arma::vec rotm2quat(const arma::mat &rotm, bool is_rotm_orthonormal)
{
  arma::vec quat(4);

  Eigen::Map<const Eigen::Matrix3d> rotm_wrapper(rotm.memptr());
  Eigen::Map<Eigen::Vector4d> quat_wrapper(quat.memptr());
  quat_wrapper = rotm2quat(rotm_wrapper, is_rotm_orthonormal);

  return quat;
}

Eigen::Vector4d mat2quat(Eigen::Matrix3d R)
{
  Eigen::Matrix4d K;
  // Calculate all elements of symmetric K matrix
  K(0,0) = R(0,0) - R(1,1) - R(2,2);
  K(0,1) = K(1,0) = R(0,1) + R(1,0);
  K(0,2) = K(2,0) = R(0,2) + R(2,0);
  K(0,3) = K(3,0) = R(2,1) - R(1,2);

  K(1,1) = R(1,1) - R(0,0) - R(2,2);
  K(1,2) = K(2,1) = R(1,2) + R(2,1);
  K(1,3) = K(3,1) = R(0,2) - R(2,0);

  K(2,2) = R(2,2) - R(0,0) - R(1,1);
  K(2,3) = K(3,2) = R(1,0) - R(0,1);

  K(3,3) = R(0,0) + R(1,1) + R(2,2);

  K = K/3;

  // For each input rotation matrix, calculate the corresponding eigenvalues
  // and eigenvectors. The eigenvector corresponding to the largest eigenvalue
  // is the unit quaternion representing the same rotation.

  Eigen::EigenSolver<Eigen::Matrix4d> es(K);

  Eigen::Vector4cd eigVal = es.eigenvalues();
  Eigen::Matrix4d eigVec = es.eigenvectors().real(); // keep only real part
  int maxIdx;
  eigVal.real().maxCoeff(&maxIdx);

  Eigen::Vector4d quat;
  quat << eigVec(3,maxIdx), eigVec(0,maxIdx), eigVec(1,maxIdx), eigVec(2,maxIdx);

  // By convention, always keep scalar quaternion element positive.
  // Note that this does not change the rotation that is represented
  // by the unit quaternion, since q and -q denote the same rotation.
  if (quat(0) < 0) quat = -quat;

  return quat;
}

// =======================================================
// =======================================================

Eigen::Matrix3d quat2rotm(Eigen::Vector4d quat)
{
  double qw=quat(0), qx=quat(1), qy=quat(2), qz=quat(3);

  Eigen::Matrix3d rotm;
  rotm << 1 - 2*qy*qy - 2*qz*qz, 	2*qx*qy - 2*qz*qw, 	2*qx*qz + 2*qy*qw,
	    2*qx*qy + 2*qz*qw, 	      1 - 2*qx*qx - 2*qz*qz, 	2*qy*qz - 2*qx*qw,
	    2*qx*qz - 2*qy*qw, 	        2*qy*qz + 2*qx*qw, 	1 - 2*qx*qx - 2*qy*qy;

  return rotm;
}

arma::mat quat2rotm(const arma::vec &quat)
{
  double qw=quat(0), qx=quat(1), qy=quat(2), qz=quat(3);

  arma::mat rotm;
  rotm = {{1 - 2*qy*qy - 2*qz*qz,      2*qx*qy - 2*qz*qw,      2*qx*qz + 2*qy*qw},
	        {    2*qx*qy + 2*qz*qw,  1 - 2*qx*qx - 2*qz*qz,      2*qy*qz - 2*qx*qw},
	        {    2*qx*qz - 2*qy*qw,      2*qy*qz + 2*qx*qw,  1 - 2*qx*qx - 2*qy*qy}};
  // rotm << 1 - 2*qy*qy - 2*qz*qz << 	2*qx*qy - 2*qz*qw     <<  	2*qx*qz + 2*qy*qw << arma::endr
	//      << 2*qx*qy + 2*qz*qw     <<  1 - 2*qx*qx - 2*qz*qz <<  	2*qy*qz - 2*qx*qw << arma::endr
	//      << 2*qx*qz - 2*qy*qw     <<    2*qy*qz + 2*qx*qw   << 	1 - 2*qx*qx - 2*qy*qy;

  return rotm;
}

// =======================================================
// =======================================================

Eigen::Vector4d rotm2axang(Eigen::Matrix3d rotm)
{
  Eigen::AngleAxis<double> angleAxis(rotm);

  Eigen::Vector4d axang;
  axang(3) = angleAxis.angle();
  axang.segment(0,3) = angleAxis.axis();

  return axang;
}

arma::vec rotm2axang(const arma::mat &rotm)
{
  arma::vec axang(4);

  Eigen::Map<const Eigen::Matrix3d> rotm_wrapper(rotm.memptr());
  Eigen::Map<Eigen::Vector4d> axang_wrapper(axang.memptr());
  axang_wrapper = rotm2axang(rotm_wrapper);

  return axang;
}

// =======================================================
// =======================================================

Eigen::Matrix3d axang2rotm(Eigen::Vector4d axang)
{
  Eigen::Matrix3d rotm;
  Eigen::Vector3d axis = axang.segment(0,3);
  axis /= axis.norm();
  double angle = axang(3);

  double x=axis(0), y=axis(1), z=axis(2), c=std::cos(angle), s=std::sin(angle), t=1-c;
  rotm <<   t*x*x + c,	    t*x*y - z*s,     t*x*z + y*s,
  	        t*x*y + z*s,    t*y*y + c,	     t*y*z - x*s,
            t*x*z - y*s,    t*y*z + x*s,     t*z*z + c;

  return rotm;
}

arma::mat axang2rotm(const arma::vec &axang)
{
  arma::vec axis = axang.subvec(0,2);
  axis /= arma::norm(axis);
  double angle = axang(3);
  double x=axis(0), y=axis(1), z=axis(2), c=std::cos(angle), s=std::sin(angle), t=1-c;

  return { { t*x*x + c,	    t*x*y - z*s,   t*x*z + y*s },
  	       { t*x*y + z*s,   t*y*y + c,	   t*y*z - x*s },
           { t*x*z - y*s,   t*y*z + x*s,   t*z*z + c   } };
}

// =======================================================
// =======================================================

Eigen::Vector4d axang2quat(Eigen::Vector4d axang)
{
  Eigen::Vector4d quat;
  double theta = axang(3);

  quat(0) = std::cos(theta/2);
  quat.segment(1,3) = std::sin(theta/2) * axang.segment(0,3);

  return quat;
}

arma::vec axang2quat(const arma::vec &axang)
{
  arma::vec quat(4);

  Eigen::Map<const Eigen::Vector4d> axang_wrapper(axang.memptr());
  Eigen::Map<Eigen::Vector4d> quat_wrapper(quat.memptr());
  quat_wrapper = axang2quat(axang_wrapper);

  return quat;
}

// =======================================================
// =======================================================

Eigen::Vector4d quat2axang(Eigen::Vector4d quat)
{
  Eigen::Vector4d axang;
  Eigen::Vector3d r = quat.segment(1,3);

  if (r.norm()){
    axang(3) = 2 * std::acos(quat(0));
    axang.segment(0,3) = r/r.norm();
  }else axang << 0, 0, 1, 0;

  return axang;
}

arma::vec quat2axang(const arma::vec &quat)
{
  arma::vec axang(4);

  Eigen::Map<const Eigen::Vector4d> quat_wrapper(quat.memptr());
  Eigen::Map<Eigen::Vector4d> axang_wrapper(axang.memptr());
  axang_wrapper = quat2axang(quat_wrapper);

  return axang;
}

// =======================================================
// =======================================================

} // namespace gmp_

} // namespace as64_
