#include <gmp_lib/GMP/GMPo.h>
#include <gmp_lib/math/quaternions.h>

namespace as64_
{

namespace gmp_
{

  const long double GMPo::zero_tol = 1e-8;

// class GMPo : public GMP

// public:

  GMPo::GMPo(unsigned N_kernels, double kernels_std_scaling): GMP(3, N_kernels, kernels_std_scaling)
  {
    this->setScaleMethod( gmp_::TrajScale::Ptr(new TrajScale_Rot_min()) );

    this->Qd0 = {1, 0, 0, 0};
    this->setQ0({1, 0, 0, 0});

    this->K = 10 * arma::vec().ones(3);
    this->D = 2 * arma::vec().ones(3);
  }

  void GMPo::train(const std::string &train_method, const arma::rowvec &Time, const arma::mat &Quat_data, arma::vec *train_error, arma::mat *Sw)
  {
    unsigned n_data = Time.size();
    this->setQ0(Quat_data.col(0));

    this->Qd0 = Quat_data.col(0);

    arma::mat qd_data(3, n_data);
    for (int j=0; j<n_data; j++) qd_data.col(j) = GMPo::quat2q(Quat_data.col(j), this->Q0);

    GMP::train(train_method, Time, qd_data, train_error, Sw);
  }

  void GMPo::setQ0(const arma::vec &Q0)
  {
    this->Q0 = Q0;
  }

  void GMPo::setQg(const arma::vec &Qg)
  {
    arma::vec qg = GMPo::quat2q(Qg, this->Q0);
    this->setGoal(qg);
  }

  arma::vec GMPo::getQd(double x) const
  {
    return GMPo::q2quat(this->getYd(x), this->Q0);
  }

  arma::vec GMPo::getVd(double x, double x_dot) const
  {
    //Ts = 0.002;
    //Qd = this->getQd(x);
    //Qd2 = this->getQd(x+x_dot*Ts);
    //vd = gmp_::quatLog(quatDiff(Qd2,Qd)) / Ts;

    arma::vec Q1 = gmp_::quatExp(this->getYd(x));
    return qLogDot_to_rotVel(this->getYdDot(x,x_dot), Q1);

  }

  arma::vec GMPo::getVdDot(double x, double x_dot, double x_ddot) const
  {
    //Ts = 0.002;
    //vd = this->getVd(x, x_dot);
    //vd2 = this->getVd(x + x_dot*Ts, x_dot + x_ddot*Ts);
    //vd_dot = (vd2 - vd) / Ts;

    arma::vec Q1 = gmp_::quatExp(this->getYd(x));
    arma::vec vd = qLogDot_to_rotVel(this->getYdDot(x,x_dot), Q1);
    return qLogDDot_to_rotAccel(this->getYdDDot(x,x_dot,x_ddot), vd, Q1);
  }

  void GMPo::getRefTraj(double x, double x_dot, double x_ddot, arma::vec &Qd, arma::vec &vd, arma::vec &vd_dot) const
  {
    arma::vec yd = this->getYd(x);
    arma::vec Q1 = gmp_::quatExp(this->getYd(x));

    Qd = GMPo::q2quat(yd, this->Q0);
    vd = qLogDot_to_rotVel(this->getYdDot(x,x_dot), Q1);
    vd_dot = qLogDDot_to_rotAccel(this->getYdDDot(x,x_dot,x_ddot), vd, Q1);
  }

  arma::vec GMPo::getRotVel(const arma::vec &Q) const
  {
    arma::vec Q1 = GMPo::getQ1(Q, this->Q0);
    return qLogDot_to_rotVel(this->getYdot(), Q1);
  }

  arma::vec GMPo::getRotAccel(const arma::vec &Q, arma::vec yc_dot) const
  {
    arma::vec Q1 = GMPo::getQ1(Q, this->Q0);
    arma::vec qddot = this->getYddot(yc_dot);
    arma::vec rotVel = qLogDot_to_rotVel(this->getYdot(), Q1);
    return qLogDDot_to_rotAccel(qddot, rotVel, Q1);
  }

  arma::vec GMPo::calcRotAccel(const gmp_::Phase &s, const arma::vec &Q, const arma::vec &rotVel, arma::vec y_c, arma::vec z_c, arma::vec yc_dot)
  {
    arma::vec Q1 = GMPo::getQ1(Q, this->Q0);
    arma::vec q = gmp_::quatLog(Q1); // q = GMPo::quat2q(Q, this->Q0);
    arma::vec invQ1 = gmp_::quatInv(Q1);

    arma::mat JQq = jacob_Q_qLog(Q1);
    arma::mat JQq_dot = jacobDot_Q_qLog(Q1, rotVel);

    arma::vec qdot = rotVel_to_qLogDot(rotVel, Q1);
    arma::vec qddot = this->calcYddot(s, q, qdot, y_c, z_c, yc_dot);

    arma::vec rotAccel = 2*gmp_::quatProd(JQq_dot*qdot + JQq*qddot, invQ1);
    return rotAccel.subvec(1,3);
  }

  arma::vec GMPo::getY(const arma::vec &Q) const
  {
    return GMPo::quat2q(Q, this->Q0);
  }

  arma::vec GMPo::getZ(const arma::vec &rotVel, const arma::vec &Q) const
  {
    arma::vec Q1 = GMPo::getQ1(Q, this->Q0);
    arma::vec qdot = rotVel_to_qLogDot(rotVel, Q1);
    arma::vec z = qdot;
    return z;
  }

  // ===========================================================
  // ***********************************************************
  // ************ Static Public Functions      ************
  // ***********************************************************
  // ===========================================================

  arma::vec GMPo::getQ1(const arma::vec &Q, const arma::vec &Q0)
  {
    return quatProd(Q, quatInv(Q0));
    // return Q;
  }


  arma::vec GMPo::quat2q(const arma::vec &Q, const arma::vec &Q0)
  {
    return quatLog(GMPo::getQ1(Q,Q0));
  }


  arma::vec GMPo::q2quat(const arma::vec &q, const arma::vec &Q0)
  {
    return quatProd( quatExp(q), Q0);
    // return quatExp(q);
  }

} // namespace gmp_

} // namespace as64_
