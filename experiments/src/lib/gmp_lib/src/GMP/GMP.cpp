#include <gmp_lib/GMP/GMP.h>
#include <gmp_lib/GMP/GMPo_Update.h>

namespace as64_
{

namespace gmp_
{

  /** GMP constructor.
   *  @param[in] n_dofs: number of degrees of freedom.
   *  @param[in] N_kernels: the number of kernels
   *  @param[in] kern_std_scale: Scaling for std of kernels (optional, default=1).
   */
  GMP::GMP(unsigned n_dofs, unsigned N_kernels, double kern_std_scale) : GMP_regressor(N_kernels, kern_std_scale)
  {
    this->K = 100 * arma::vec().ones(n_dofs);
    this->D = 30 * arma::vec().ones(n_dofs);

    this->W = arma::mat().zeros(n_dofs, N_kernels);

    this->Y0d = arma::vec().zeros(n_dofs);
    this->Ygd = arma::vec().ones(n_dofs);
    this->Y0 = this->Y0d;
    this->Yg = this->Ygd;

    this->setScaleMethod( gmp_::TrajScale::Ptr(new TrajScale_Prop(n_dofs)) );

    this->y_dot = arma::vec().zeros(n_dofs);
    this->z_dot = arma::vec().zeros(n_dofs);
  }

  /** Returns the number of DoFs.
   *  return: number of DoFs.
   */
  unsigned GMP::numOfDoFs() const
  {
    return this->W.n_rows;
  }

  /** Returns the number of kernels.
   *  @return: number of kernels.
   */
  unsigned GMP::numOfKernels() const
  {
    return this->c.size();
  }

  /** Sets the initial position.
   *  @param[in] y0: initial position.
   */
  void GMP::setY0(const arma::vec &y0)
  {
    this->Y0 = y0;
    this->traj_sc->setNewStartFinalPos(this->Y0, this->Yg);
  }

  arma::vec GMP::getY0() const
  {
    return this->Y0;
  }

  arma::vec GMP::getY0d()
  {
    return this->Y0d;
  }


  /** Set goal position.
   *  @param[in] g: goal position.
   */
  void GMP::setGoal(const arma::vec &g)
  {
    this->Yg = g;
    this->traj_sc->setNewStartFinalPos(this->Y0, this->Yg);
  }

  arma::vec GMP::getGoal() const
  {
    return this->Yg;
  }

  arma::mat GMP::getScaling() const
  {
    return this->traj_sc->getScaling();
  }

  arma::mat GMP::getInvScaling() const
  {
    return this->traj_sc->getInvScaling();
  }

  void GMP::setScaleMethod(const std::shared_ptr<gmp_::TrajScale> &traj_scale_obj)
  {
      if (this->numOfDoFs() != traj_scale_obj->numOfDoFs())
        throw std::runtime_error("\33[1m\33[31m[GMP::setScaleMethod]: Incompatible number of DoFs...\33[0m");

      this->traj_sc = traj_scale_obj;
      this->traj_sc->setNominalStartFinalPos(this->Y0d, this->Ygd);
      this->traj_sc->setNewStartFinalPos(this->Y0, this->Yg);
  }

  void GMP::train(const std::string &train_method, const arma::rowvec &x, const arma::mat &yd_data, arma::vec *train_err, arma::mat *Sw)
  {
    for (int i=0; i<x.size(); i++)
    { if (x(i)>1 || x(i)<0) throw std::runtime_error("\33[1m\33[31m[GMP::train]: The training timestamps are not normalized...\33[0m"); }

    unsigned n_data = x.size();
    unsigned num_ker = this->numOfKernels();
    unsigned n_dofs = this->numOfDoFs();
    unsigned i_end = n_data-1;

    arma::mat H(num_ker, n_data);
    for (int j=0; j<n_data; j++) H.col(j) = this->regressVec(x(j));

    if (train_method.compare("LWR") == 0)
      this->W = yd_data*H.t() / arma::repmat(arma::sum(H,1).t(),n_dofs,1);
    else if (train_method.compare("LS") == 0)
      this->W = arma::solve(H.t(), yd_data.t()).t(); // yd_data / H;
    else
      throw std::runtime_error("\33[1m\33[31m[GMP::train]: Unsupported training method...\33[0m");

    this->Y0d = this->W*H.col(0);
    this->Ygd = this->W*H.col(i_end);

    this->setY0(this->Y0d);
    this->setGoal(this->Ygd);

    this->traj_sc->setNominalStartFinalPos(this->Y0d, this->Ygd);
    this->traj_sc->setNewStartFinalPos(this->getY0(), this->getGoal());

    if (train_err)
    {
        *train_err = arma::vec().zeros(n_dofs);
        arma::mat err_data = this->W*H - yd_data;
        for (int i=0; i<n_dofs; i++) train_err->at(i) = arma::norm(err_data.row(i));
    }

    if (Sw) *Sw = arma::inv(H*H.t());

  }

  void GMP::autoRetrain(unsigned N_kernels, double kern_std_scale, unsigned n_points, const std::string &train_method, arma::vec *train_err, arma::mat *Sw)
  {     
      arma::rowvec x = arma::linspace<arma::rowvec>(0,1, n_points);
      arma::mat yd_data(this->numOfDoFs(), x.size());
      for (int j=0;j<x.size();j++) yd_data.col(j) = this->W*this->regressVec(x(j));
      
      this->setKernels(N_kernels, kern_std_scale);
      this->train(train_method, x, yd_data, train_err, Sw);
  }

  void GMP::update(const gmp_::Phase &s, const arma::vec &y, const arma::vec &z, arma::vec y_c, arma::vec z_c)
  {
    unsigned n_dofs = this->numOfDoFs();

    if (y_c.size() == 1) y_c = arma::vec().ones(n_dofs)*y_c(0);
    if (z_c.size() == 1) z_c = arma::vec().ones(n_dofs)*z_c(0);

    arma::vec yd = this->getYd(s.x);
    arma::vec yd_dot = this->getYdDot(s.x, s.x_dot);
    arma::vec yd_ddot = this->getYdDDot(s.x, s.x_dot, s.x_ddot);

    this->y_dot = z + y_c;
    this->z_dot = this->K%(yd - y) + this->D%(yd_dot - z) + yd_ddot + z_c;
  }

  arma::vec GMP::getYdot() const
  {
    return this->y_dot;
  }

  arma::vec GMP::getZdot() const
  {
    return this->z_dot;
  }

  arma::vec GMP::getYddot(arma::vec yc_dot) const
  {
    unsigned n_dofs = this->numOfDoFs();
    if (yc_dot.size()==1) yc_dot = arma::vec().ones(n_dofs)*yc_dot(0);

    return this->getZdot() + yc_dot;
  }

  arma::vec GMP::calcYddot(const gmp_::Phase &s, const arma::vec &y, const arma::vec &y_dot, arma::vec y_c, arma::vec z_c, arma::vec yc_dot) const
  {
      unsigned n_dofs = this->numOfDoFs();

      if (y_c.size()==1) y_c = arma::vec().ones(n_dofs)*y_c(0);
      if (z_c.size()==1) z_c = arma::vec().ones(n_dofs)*z_c(0);
      if (yc_dot.size()==1) yc_dot = arma::vec().ones(n_dofs)*yc_dot(0);

      arma::vec yd = this->getYd(s.x);
      arma::vec yd_dot = this->getYdDot(s.x, s.x_dot);
      arma::vec yd_ddot = this->getYdDDot(s.x, s.x_dot, s.x_ddot);

      arma::vec z = y_dot - y_c;
      arma::vec z_dot = this->K%(yd - y) + this->D%(yd_dot - z) + yd_ddot + z_c;
      return (z_dot + yc_dot);
  }

  arma::vec GMP::getYd(double x) const
  {
    return this->getScaling()*(this->W*this->regressVec(x) - this->Y0d) + this->Y0;
  }

  arma::vec GMP::getYdDot(double x, double x_dot) const
  {
    return this->getScaling()*this->W*this->regressVecDot(x,x_dot);
  }


  arma::vec GMP::getYdDDot(double x, double x_dot, double x_ddot) const
  {
    return this->getScaling()*this->W*this->regressVecDDot(x,x_dot,x_ddot);
  }

  void GMP::deepCopy(gmp_::GMP *cp_obj) const
  {
    // make a shallow copy first
    *cp_obj = *this;

    // make a deep copy of pointer-like objects
    cp_obj->traj_sc = this->traj_sc->deepCopy();
  }

} // namespace gmp_

} // namespace as64_
