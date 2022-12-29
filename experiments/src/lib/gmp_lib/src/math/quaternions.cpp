#include <gmp_lib/math/quaternions.h>
#include <cmath>

namespace as64_
{

namespace gmp_
{

Eigen::Matrix3d vec2ssMat(const Eigen::Vector3d &v)
{
  Eigen::Matrix3d ssMat;
  ssMat(0, 1) = -v(2);
  ssMat(0, 2) =  v(1);
  ssMat(1, 0) =  v(2);
  ssMat(1, 2) = -v(0);
  ssMat(2, 0) = -v(1);
  ssMat(2, 1) =  v(0);

  return ssMat;
}

arma::mat vec2ssMat(const arma::vec &v)
{
  arma::mat ssMat(3,3);
  ssMat(0, 1) = -v(2);
  ssMat(0, 2) =  v(1);
  ssMat(1, 0) =  v(2);
  ssMat(1, 2) = -v(0);
  ssMat(2, 0) = -v(1);
  ssMat(2, 1) =  v(0);

  return ssMat;
}

Eigen::Vector4d quatInv(const Eigen::Vector4d &quat)
{
  Eigen::Vector4d quatI;

  quatI(0) = quat(0);
  quatI.segment(1,3) = - quat.segment(1,3);

  return quatI;
}

arma::vec quatInv(const arma::vec &quat)
{
  arma::vec quatI(4);

  quatI(0) = quat(0);
  quatI.subvec(1,3) = - quat.subvec(1,3);

  return quatI;
}

Eigen::Matrix4d quat2mat(const Eigen::Vector4d &quat)
{
  Eigen::Matrix4d mat;
  Eigen::Vector3d e = quat.segment(1,3);
  double n = quat(0);

  mat << n, -e.transpose(), e, n*Eigen::Matrix3d::Identity() + vec2ssMat(e);

  return mat;
}

Eigen::Vector4d quatProd(const Eigen::Vector4d &quat1, const Eigen::Vector4d &quat2)
{
  Eigen::Vector4d quat12;

  // quat12 = quat2qmat(quat1) * quat2;
  double n1 = quat1(0);
  Eigen::Vector3d e1 = quat1.segment(1,3);

  double n2 = quat2(0);
  Eigen::Vector3d e2 = quat2.segment(1,3);

  quat12(0) = n1*n2 - e1.dot(e2);
  quat12.segment(1,3) = n1*e2 + n2*e1 + e1.cross(e2);

  return quat12;
}

arma::vec quatProd(const arma::vec &quat1, const arma::vec &quat2)
{
  arma::vec quat12(4);

  double n1 = quat1(0);
  arma::vec e1 = quat1.subvec(1,3);

  double n2 = quat2(0);
  arma::vec e2 = quat2.subvec(1,3);

  quat12(0) = n1*n2 - arma::dot(e1,e2);
  quat12.subvec(1,3) = n1*e2 + n2*e1 + arma::cross(e1,e2);

  return quat12;
}

Eigen::Vector4d quatExp(const Eigen::Vector3d &v_rot, double zero_tol)
{
  Eigen::Vector4d quat;
  double norm_v_rot = v_rot.norm();
  double theta = norm_v_rot;

 if (norm_v_rot > zero_tol)
 {
    quat(0) = std::cos(theta/2);
    quat.segment(1,3) = std::sin(theta/2)*v_rot/norm_v_rot;
  }
  else{
    quat << 1, 0, 0, 0;
  }

  return quat;
}

arma::vec quatExp(const arma::vec &v_rot, double zero_tol)
{
  arma::vec quat(4);
  double norm_v_rot = arma::norm(v_rot);
  double theta = norm_v_rot;

 if (norm_v_rot > zero_tol)
 {
    quat(0) = std::cos(theta/2);
    quat.subvec(1,3) = std::sin(theta/2)*v_rot/norm_v_rot;
  }
  else{
    quat << 1 << 0 << 0 << 0;
  }

  return quat;
}

Eigen::Vector3d quatLog(const Eigen::Vector4d &quat, double zero_tol)
{
  Eigen::Vector3d e = quat.segment(1,3);
  double n = quat(0);

  Eigen::Vector3d q_log;
  double norm_e = e.norm();

  if (norm_e > zero_tol) q_log = 2*std::atan2(norm_e,n)*e/norm_e;
  else q_log = Eigen::Vector3d::Zero();

  return q_log;
}

arma::vec quatLog(const arma::vec &quat, double zero_tol)
{
  arma::vec e = quat.subvec(1,3);
  double n = quat(0);

  if (n > 1) n = 1;
  if (n < -1) n = -1;

  arma::vec q_log(3);
  double norm_e = arma::norm(e);

  if (norm_e > zero_tol) q_log = 2*std::atan2(norm_e,n)*e/norm_e;
  else q_log = arma::vec().zeros(3);

  return q_log;
}

arma::vec quatDiff(const arma::vec &quat1, const arma::vec &quat2)
{
  arma::vec quatD;

  arma::vec Q2 = quat2;
  if (arma::dot(quat1,Q2)<0) Q2 = -Q2;

  quatD = quatProd(quat1, quatInv(Q2));

  return quatD;
}

arma::vec quatLogDiff(const arma::vec &quat1, const arma::vec &quat2)
{
  return gmp_::quatLog(gmp_::quatDiff(quat1,quat2));
}

// ===========================================================
// ===========================================================

arma::vec rotVel_to_qLogDot(const arma::vec &rotVel, const arma::vec &logQ)
{
  arma::mat JqQ = jacob_qLog_Q(logQ);
  return 0.5 * JqQ * quatProd( arma::join_vert(arma::vec({0}), rotVel), logQ );
}


arma::vec qLogDot_to_rotVel(const arma::vec &logQ_dot, const arma::vec &Q)
{
  arma::mat JQq = jacob_Q_qLog(Q);
  arma::vec rotVel = 2 * quatProd( JQq*logQ_dot, quatInv(Q) );
  return rotVel.subvec(1,3);
}


arma::vec rotAccel_to_qLogDDot(const arma::vec &rotAccel, const arma::vec &rotVel, const arma::vec &Q)
{
  arma::vec rotVelQ = arma::join_vert(arma::vec({0}), rotVel);
  arma::vec rotAccelQ = arma::join_vert(arma::vec({0}), rotAccel);

  arma::mat J = jacob_qLog_Q(Q);
  arma::mat Jdot = jacobDot_qLog_Q(Q, rotVel);

  return 0.5*(Jdot * quatProd(rotVelQ, Q) + J * quatProd( rotAccelQ+0.5*quatProd(rotVelQ,rotVelQ), Q ) );
}


arma::vec qLogDDot_to_rotAccel(const arma::vec &logQ_ddot, const arma::vec &rotVel, const arma::vec &Q)
{
  arma::vec logQ_dot = rotVel_to_qLogDot(rotVel, Q);
  arma::mat Jq = jacob_Q_qLog(Q);
  arma::mat dJq = jacobDot_Q_qLog(Q, rotVel);

  arma::vec rotAccel = 2 * ( quatProd( dJq*logQ_dot + Jq*logQ_ddot, quatInv(Q) ) );
  return rotAccel.subvec(1,3);
}

const long double JLOG_ZERO_TOL = 1e-8;

arma::mat jacob_Q_qLog(const arma::vec &Q)
{
  arma::mat JQq(4,3);

  if ( (1-std::fabs(Q(0))) <= JLOG_ZERO_TOL)
  {
    JQq.row(0) = arma::rowvec().zeros(3);
    JQq.submat(1,0,3,2) = arma::mat().eye(3,3);
    return JQq;
  }

  double w = Q(0);
  arma::vec v = Q.subvec(1,3);
  double norm_v = arma::norm(v);
  arma::vec eta = v / norm_v;
  double s_th = norm_v;
  double c_th = w;
  double th = std::atan2(s_th, c_th);
  arma::mat Eta = eta*eta.t();

  JQq.row(0) = -0.5 * s_th * eta.t();
  JQq.submat(1,0,3,2) = 0.5 * ( (arma::mat().eye(3,3) - Eta)*s_th/th + c_th*Eta );
  return JQq;
}

arma::mat jacob_qLog_Q(const arma::vec &Q)
{
  arma::mat JqQ(3,4);

  if ( (1-std::fabs(Q(0))) <= JLOG_ZERO_TOL)
  {
    JqQ.col(0) = arma::vec().zeros(3);
    JqQ.submat(0,1,2,3) = arma::mat().eye(3,3);
    return JqQ;
  }

  double w = Q(0);
  arma::vec v = Q.subvec(1,3);
  double norm_v = arma::norm(v);
  arma::vec eta = v / norm_v;
  double s_th = norm_v;
  double c_th = w;
  double th = std::atan2(s_th, c_th);

  JqQ.col(0) = 2*eta*(th*c_th - s_th)/std::pow(s_th,2);
  JqQ.submat(0,1,2,3) = 2*arma::mat().eye(3,3)*th/s_th;
  return JqQ;
}


arma::mat jacobDot_qLog_Q(const arma::vec &Q, const arma::vec &rotVel)
{
  arma::mat JqQ_dot(3,4);

  arma::vec qdot = rotVel_to_qLogDot(rotVel, Q);

  if ( (1-std::fabs(Q(0))) <= JLOG_ZERO_TOL)
  {
    JqQ_dot.col(0) = -qdot/3;
    JqQ_dot.submat(0,1,2,3) = arma::mat().zeros(3,3);
    return JqQ_dot;
  }

  double w = Q(0);
  arma::vec v = Q.subvec(1,3);
  double norm_v = arma::norm(v);
  arma::vec eta = v / norm_v;
  double s_th = norm_v;
  double c_th = w;
  double th = std::atan2(s_th, c_th);
  arma::mat Eta = eta*eta.t();
  double temp = (th*c_th-s_th)/std::pow(s_th,2);

  JqQ_dot.col(0) = ((-th/s_th - 2*c_th*temp/s_th)*Eta + temp*(arma::mat().eye(3,3)-Eta)/th)*qdot;
  JqQ_dot.submat(0,1,2,3) = (-temp*arma::dot(eta,qdot))*arma::mat().eye(3,3);
  return JqQ_dot;
}


arma::mat jacobDot_Q_qLog(const arma::vec &Q, const arma::vec &rotVel)
{
  arma::mat JQq_dot(4,3);

  arma::vec qdot = rotVel_to_qLogDot(rotVel, Q);

  if ( (1-std::fabs(Q(0))) <= JLOG_ZERO_TOL)
  {
    JQq_dot.row(0) = -qdot.t()/4;
    JQq_dot.submat(1,0,3,2) = arma::mat().zeros(3,3);
    return JQq_dot;
  }

  double w = Q(0);
  arma::vec v = Q.subvec(1,3);
  double norm_v = arma::norm(v);
  arma::vec eta = v / norm_v;
  double s_th = norm_v;
  double c_th = w;
  double th = std::atan2(s_th, c_th);
  arma::mat Eta = eta*eta.t();
  arma::mat I_eta = arma::mat().eye(3,3) - Eta;
  double temp = ((th*c_th-s_th)/std::pow(th,2));

  JQq_dot.row(0) = -0.25 * qdot.t() * (c_th*Eta + (s_th/th)*I_eta);
  JQq_dot.submat(1,0,3,2) = (0.25*arma::dot(eta,qdot))*( temp*I_eta - s_th*Eta ) + 0.25*temp*( eta*(qdot.t()*I_eta) + (I_eta*qdot)*eta.t() );

  return JQq_dot;
}


} // namespace gmp_

} // namespace as64_
