#include <online_adapt_controller/utils/DMP_pp.h>

// #include <gmp_lib/gmp_lib.h>
#include <gmp_lib/io/gmp_io.h>
#include <gmp_lib/math/math.h>

// #include <ros/package.h>

using namespace as64_;

#define DMP_pp_fun_ std::string("[online_adapt_::DMP_pp::") + __func__ + "]: "

DMP_pp::DMP_pp()
{
  initialized = false;
}

DMP_pp::~DMP_pp()
{

}

void DMP_pp::reset(bool reset_model)
{
  // if (first_exec || reload_model) this->load(gmp_model_filename);
  // first_exec = false;
  if (!initialized) throw std::runtime_error(DMP_pp_fun_ + "The model is not initialized. Call \"load()\" first!");
  
  if (reset_model)
  {
    gmp_p->W = W0_p;
    gmp_o->W = W0_o;
  }

  vp_map.clear();

  gmp_p_up.reset(new gmp_::GMP_Update(gmp_p.get()));
  gmp_p_up->syncUpdate(sync_update);
  gmp_p_up->recursiveUpdate(recursive_update);

  gmp_o_up.reset(new gmp_::GMPo_Update(gmp_o.get()));
  gmp_o_up->syncUpdate(sync_update);
  gmp_o_up->recursiveUpdate(recursive_update);

  resetSigmaw();

  prev_pos_state.initialized = false;
  prev_quat_state.initialized = false;
}

void DMP_pp::resetSigmaw()
{
  arma::rowvec s_data = arma::linspace<arma::rowvec>(0, 1, 100);
  if (Sigma_w_type.compare("pos") == 0)
  {
    gmp_p_up->initSigmaWfromPosMsr(s_data, Sigma0_tol);
    gmp_o_up->initSigmaWfromPosMsr(s_data, Sigma0_tol);
  }
  else if (Sigma_w_type.compare("vel") == 0)
  {
    gmp_p_up->initSigmaWfromVelMsr(s_data, 1.0, Sigma0_tol);
    gmp_o_up->initSigmaWfromVelMsr(s_data, 1.0, Sigma0_tol);
  }
  else if (Sigma_w_type.compare("accel") == 0)
  {
    gmp_p_up->initSigmaWfromAccelMsr(s_data, 1.0, Sigma0_tol);
    gmp_o_up->initSigmaWfromAccelMsr(s_data, 1.0, Sigma0_tol);
  }
  else if (Sigma_w_type.compare("exp") == 0)
  {
    gmp_p_up->initExpSigmaw(0.01);
    gmp_o_up->initExpSigmaw(0.01);
  }
  else throw std::runtime_error("Unrecognized \"Sigma_w_type\"...");
}

// ============= Update Model ================

void DMP_pp::initUpdate(double s_start, const arma::vec &P0, const arma::vec &Pg, const arma::vec &Q0, const arma::vec &Qg, double Tf)
{
  double s0 = 0;
  double sf = 1;
  double s_dot = 1 / Tf;
  double s_ddot = 0;

  arma::vec O_zeros = arma::vec().zeros(3);

  gmp_p->setY0(P0);
  gmp_p->setGoal(Pg);
  gmp_p->setScaleMethod(std::shared_ptr<gmp_::TrajScale>(new gmp_::TrajScale_None(3)));

  gmp_o->setQ0(Q0);
  gmp_o->setQg(Qg);
  gmp_o->setScaleMethod(std::shared_ptr<gmp_::TrajScale>(new gmp_::TrajScale_None(3)));

  // ---------- Position ------------

  // initial state constraints: pos, vel, accel
  gmp_p_up->updatePos(s0, P0, r0(0));
  gmp_p_up->updateVel(s0, s_dot, O_zeros, r0(1));
  gmp_p_up->updateAccel(s0, s_dot, s_ddot, O_zeros, r0(2));

  // final state constraints: pos, vel, accel
  gmp_p_up->updatePos(sf, Pg, rf(0));
  gmp_p_up->updateVel(sf, s_dot, O_zeros, rf(1));
  gmp_p_up->updateAccel(sf, s_dot, s_ddot, O_zeros, rf(2));

  gmp_p_up->updateNow();

  prev_pos_state.P0 = P0;
  prev_pos_state.Pg = Pg;
  prev_pos_state.s = s_start;
  prev_pos_state.s_dot = s_dot;
  prev_pos_state.initialized = true;

  // ---------- Orientation ------------

  // initial state constraints: pos, vel, accel
  gmp_o_up->updateQuat(s0, Q0, r0(0));
  gmp_o_up->updateRotVel(s0, s_dot, O_zeros, Q0, r0(1));
  gmp_o_up->updateRotAccel(s0, s_dot, s_ddot, O_zeros, O_zeros, Q0, r0(2));

  // final state constraints: pos, vel, accel
  gmp_o_up->updateQuat(sf, Qg, rf(0));
  gmp_o_up->updateRotVel(sf, s_dot, O_zeros, Qg, rf(1));
  gmp_o_up->updateRotAccel(sf, s_dot, s_ddot, O_zeros, O_zeros, Qg, rf(2));

  gmp_o_up->updateNow();

  prev_quat_state.Q0 = Q0;
  prev_quat_state.Qg = Qg;
  prev_quat_state.s = s_start;
  prev_quat_state.s_dot = s_dot;
  prev_quat_state.initialized = true;
}

void DMP_pp::updatePos(const arma::vec &Pg, double s, double s_dot, const arma::vec &P, const arma::vec &P_dot)
{
  if (!prev_pos_state.initialized)
    throw std::runtime_error(DMP_pp_fun_ + "You need to call DMP_pp::initUpdate first!");

  double s0 = 0;
  double sf = 1;
  double s_ddot = 0;

  arma::vec P1 = P;
  arma::vec P1_dot = P_dot;
  arma::vec P1_ddot = gmp_p->getYdDDot(prev_pos_state.s, prev_pos_state.s_dot, s_ddot);

  arma::vec Pd = gmp_p->getYd(s);
  arma::vec Pd_dot = gmp_p->getYdDot(s, s_dot);

  if (!adapt_to_robot_state)
  {
    P1 = Pd;
    P1_dot = Pd_dot;
  }
  else
  {
    // Add a small tolerance to avoid numerical integration drift
    if (arma::norm(P1 - Pd) < 5e-4) P1 = Pd;
    if (arma::norm(P1_dot - Pd_dot) < 1e-3) P1_dot = Pd_dot;
  }

  // std::cerr << "ep     = " << arma::norm(P - Pd) << "\n";
  // std::cerr << "ep_dot = " << arma::norm(P1_dot - Pd_dot) << "\n";
  // std::cerr << "------------------------\n";

  arma::vec O_zeros = arma::vec().zeros(3);
  
  // ---------- downdate --------------
  if (recursive_update)  
  {
    gmp_p_up->updatePos(sf, prev_pos_state.Pg, -rf(0));
    // ===> For zero vel/accel and assuming sg_ddot=0, this is redundant:
    // gmp_p_up->updateVel(sf, prev_pos_state.s_dot, O_zeros, -rf(1));
    // gmp_p_up->updateAccel(sf, prev_pos_state.s_dot, s_ddot, O_zeros, -rf(2));

    gmp_p_up->updateNow();
  }

  // ---------- update --------------

  // current state constraints: pos, vel, accel
  gmp_p_up->updatePos(s, P1, r1(0));
  gmp_p_up->updateVel(s, s_dot, P1_dot, r1(1));
  if (enforce_cont_accel)
  {
    // gmp_p_up->updateAccel(s, s_dot, s_ddot, P1_ddot, r1(2));
    gmp_p_up->updateAccel(prev_pos_state.s, prev_pos_state.s_dot, s_ddot, P1_ddot, r1(2));
  }

  // final state constraints: pos, vel, accel
  gmp_p_up->updatePos(sf, Pg, rf(0));
  // ===> For zero vel/accel and assuming sg_ddot=0, this is redundant:
  // gmp_p_up->updateVel(sf, s_dot, O_zeros, rf(1));
  // gmp_p_up->updateAccel(sf, s_dot, s_ddot, O_zeros, rf(2));

  gmp_p_up->updateNow();

  // prev_pos_state.P0 = P0;
  prev_pos_state.Pg = Pg;
  prev_pos_state.s = s;
  prev_pos_state.s_dot = s_dot;
  prev_pos_state.initialized = true;
}

void DMP_pp::updateOrient(const arma::vec &Qg, double s, double s_dot, const arma::vec &Q, const arma::vec &vRot)
{
  if (!prev_quat_state.initialized)
    throw std::runtime_error(DMP_pp_fun_ + "You need to call DMP_pp::initUpdate first!");

  double s0 = 0;
  double sf = 1;
  double s_ddot = 0;

  arma::vec Qf = Qg;
  arma::vec Q1 = Q;
  arma::vec vRot1 = vRot;
  arma::vec q1_ddot = gmp_o->getYdDDot(prev_quat_state.s, prev_quat_state.s_dot, s_ddot);

  arma::vec Qd = gmp_o->getQd(s);
  arma::vec vRotd = gmp_o->getVd(s, s_dot);

  if (arma::dot(Q1, Qd) < 0) Q1 = -Q1;

  if (!adapt_to_robot_state)
  {
    Q1 = Qd;
    vRot1 = vRotd;
  }
  else
  {
    if (arma::norm(gmp_::quatLog(gmp_::quatDiff(Q1, Qd))) < 5e-4) Q1 = Qd;
    if (arma::norm(vRot1 - vRotd) < 1e-3) vRot1 = vRotd;
  }

  // std::cerr << "eo     = " << arma::norm(gmp_::quatLog(gmp_::quatDiff(Q1, Qd))) << "\n";
  // std::cerr << "eo_dot = " << arma::norm(vRot1 - vRotd) << "\n";
  // std::cerr << "------------------------\n";

  arma::vec O_zeros = arma::vec().zeros(3);

  if (recursive_update)  
  {
    gmp_o_up->updateQuat(sf, prev_quat_state.Qg, -rf(0));
    // gmp_o_up->updateRotVel(sf, prev_quat_state.s_dot, O_zeros, prev_quat_state.Qg, -rf(1));
    // gmp_o_up->updateRotAccel(sf, prev_quat_state.s_dot, s_ddot, O_zeros, O_zeros, prev_quat_state.Qg, -rf(2));

    gmp_o_up->updateNow();
  }

  // current state constraints: pos, vel, accel
  gmp_o_up->updateQuat(s, Q1, r1(0));
  gmp_o_up->updateRotVel(s, s_dot, vRot1, Q1, r1(1));
  if (enforce_cont_accel) gmp_o_up->updateAccel(prev_quat_state.s, prev_quat_state.s_dot, s_ddot, q1_ddot, r1(2));

  // final state constraints: pos, vel, accel
  gmp_o_up->updateQuat(sf, Qf, rf(0));
  // gmp_o_up->updateRotVel(sf, s_dot, O_zeros, Qf, rf(1));
  // gmp_o_up->updateRotAccel(sf, s_dot, s_ddot, O_zeros, O_zeros, Qf, rf(2));

  gmp_o_up->updateNow();
  
  // prev_quat_state.Q0 = Q0;
  prev_quat_state.Qg = Qf;
  prev_quat_state.s = s;
  prev_quat_state.s_dot = s_dot;
  prev_quat_state.initialized = true;
}

double DMP_pp::updateViapoint(double s_v, const arma::vec &pos, const arma::vec &quat, bool reverse)
{
  double s;
  arma::vec pos2 = pos;
  arma::vec quat2 = quat;
  
  if (reverse)
  {
    s = findClosest(s_v, 0, 60, pos);
    double ds = s*0.1 + 1e-9;
    s = findClosest(std::min(s+ds, 1.0), std::max(s_v, s-ds), 20, pos, &pos2);
  }
  else
  {
    s = findClosest(s_v, 1.0, 60, pos);
    double ds = (1-s)*0.1 + 1e-9;
    s = findClosest(std::max(s_v, s-ds), std::min(s+ds, 1.0), 20, pos, &pos2);
  }

  if (quat2.has_nan()) quat2 = this->getRefQuat(s);

  gmp_p_up->updatePos(s, pos2, rv);
  gmp_p_up->updateNow();

  gmp_o_up->updateQuat(s, quat2, rv);
  gmp_o_up->updateNow();

  return s;
}

std::vector<double> DMP_pp::updateViapoints(double s, const std::vector<arma::vec> &Pv, const std::vector<arma::vec> &Qv, const std::string &vp_set, std::vector<arma::vec> *Pv2, std::vector<arma::vec> *Qv2)
{
  if (Pv.empty()) return std::vector<double>();
  if (Pv.size() != Qv.size()) throw std::runtime_error(DMP_pp_fun_ + "Pv.size() != QV.size() (" + std::to_string(Pv.size()) + " != " + std::to_string(Qv.size()) + ")");

  std::vector<double> s_v(Pv.size());
  std::vector<arma::vec> pos_vp(Pv.size());
  std::vector<arma::vec> quat_vp(Pv.size());

  // use the last updated reference trajectory to specify the via-points' phase
  for (int i=0; i<s_v.size(); i++)
  {
    double s0 = s;
    auto &pos = Pv[i];
    arma::vec pos2 = pos;
    int n_search = std::max(int((1-s)*50+0.5), 10);
    s = findClosest(s0, 1.0, n_search, pos);
    double ds = (1-s)*0.1 + 1e-9;
    s = findClosest(std::max(s0, s-ds), std::min(s+ds, 1.0), 10, pos, &pos2);
    s_v[i] = s;

    arma::vec quat2 = Qv[i];
    if (quat2.has_nan()) quat2 = getRefQuat(s_v[i]);

    pos_vp[i] = pos2;
    quat_vp[i] = quat2;
  }

  // now downdate
  this->downdateViapoints(vp_set);

  vp_map[vp_set] = ViapointSet(s_v, pos_vp, quat_vp);

  for (int i=0; i<s_v.size(); i++)
  {
    gmp_p_up->updatePos(s_v[i], pos_vp[i], rv);
    gmp_o_up->updateQuat(s_v[i], quat_vp[i], rv);
  }

  gmp_p_up->updateNow();
  gmp_o_up->updateNow();

  // std::cerr << "sv = " << arma::rowvec(s_v) << "\n";

  if (Pv2) *Pv2 = pos_vp;
  if (Qv2) *Qv2 = quat_vp;

  return s_v;
}
  
void DMP_pp::downdateViapoints(const std::string &vp_set)
{
  auto it = vp_map.find(vp_set);
  if (it == vp_map.end())
  {
    // std::cerr << ("\33[1;33m" + DMP_pp_fun_ + "The viapoint set \"" + vp_set + "\" does not exist...\33[0m");
    return;
  }

  auto &vp = it->second;
  for (int i=0; i<vp.s.size(); i++)
  {
    gmp_p_up->updatePos(vp.s[i], vp.P[i], -rv);
    gmp_o_up->updateQuat(vp.s[i], vp.Q[i], -rv);
  }

  gmp_p_up->updateNow();
  gmp_o_up->updateNow();

  vp_map.erase(it);
}

double DMP_pp::findClosest(double s0, double sf, int n_points, const arma::vec &pos, arma::vec *p2) const
{
  arma::rowvec s_data;
  
  // choose whether to search forward or backwards in time
  if (s0 < sf) s_data = arma::linspace<arma::rowvec>(s0, sf, n_points);
  else s_data = -arma::linspace<arma::rowvec>(-s0, -sf, n_points);

  double min_dist = 1e4;
  double i_min = -1;
  for (int i=0; i<s_data.size(); i++)
  {
    double s = s_data(i);
    auto idx = arma:: find_finite(pos);
    double dist = arma::norm(pos.elem(idx) - this->getRefPos(s).elem(idx));
    if (dist < min_dist)
    {
      i_min = i;
      min_dist = dist;
    }
  }

  double s_min = s_data(i_min);

  if (p2)
  {
    *p2 = pos;
    auto idx = arma::find_nonfinite(*p2);
    p2->elem(idx) = this->getRefPos(s_min).elem(idx);
  }

  return s_min;
}


// ============= Refence Position ================

arma::vec DMP_pp::getRefPos(double s) const
{
  if (s > 1) s = 1;
  if (s < 0) s = 0;
  return gmp_p->getYd(s);
}

arma::vec DMP_pp::getRefVel(double s, double s_dot) const
{ 
  if (s > 1 || s < 0) return arma::vec().zeros(3);
  return gmp_p->getYdDot(s, s_dot);
}

arma::vec DMP_pp::getRefAccel(double s, double s_dot, double s_ddot) const
{
  if (s > 1 || s < 0) return arma::vec().zeros(3);
  return gmp_p->getYdDDot(s, s_dot, s_ddot);
}

arma::vec DMP_pp::getMotionDirection(double s) const
{
  arma::vec v = getRefVel(s, 1);
  double v_norm = arma::norm(v);
  if (v_norm < 1e-3) v = arma::vec().zeros(3);
  else v = v / v_norm;
  return v;
}
// ============= Refence Orientation ================

arma::vec DMP_pp::getRefQuat(double s) const
{
  if (s > 1) s = 1;
  if (s < 0) s = 0;
  return gmp_o->getQd(s);
}

arma::vec DMP_pp::getRefRotVel(double s, double s_dot) const
{
  if (s > 1 || s < 0) return arma::vec().zeros(3);
  return gmp_o->getVd(s, s_dot);
}

arma::vec DMP_pp::getRefRotAccel(double s, double s_dot, double s_ddot) const
{
  if (s > 1 || s < 0) return arma::vec().zeros(3);
  return gmp_o->getVdDot(s, s_dot, s_ddot);
}

// ============ Save/Load model ==============

bool DMP_pp::load(const std::string &filename, std::string *err_msg)
{
  try
  {
    gmp_::FileIO fid(filename, gmp_::FileIO::in);

    gmp_p.reset(new gmp_::GMP());
    gmp_::read(gmp_p.get(), filename, "pos_");
    gmp_p->setScaleMethod(gmp_::TrajScale::Ptr(new gmp_::TrajScale_None(gmp_p->numOfDoFs())));

    gmp_o.reset(new gmp_::GMPo());
    gmp_::read(gmp_o.get(), filename, "orient_");
    gmp_o->setScaleMethod(gmp_::TrajScale::Ptr(new gmp_::TrajScale_None(gmp_o->numOfDoFs())));

    this->W0_p = gmp_p->W;
    this->W0_o = gmp_o->W;

    initialized = true;
    
    fid.close();
    return true;
  }
  catch(const std::exception& e)
  {
    if (err_msg) *err_msg = DMP_pp_fun_ + e.what();
    return false;
  }
}

bool DMP_pp::save(const std::string &filename, std::string *err_msg)
{
  try
  {
    if (!initialized) throw std::runtime_error("The model is empty...\n");

    gmp_::FileIO fid(filename, gmp_::FileIO::out | gmp_::FileIO::trunc);
    gmp_::write(gmp_p.get(), fid, "pos_");
    gmp_::write(gmp_o.get(), fid, "orient_");
    fid.close();
    return true;
  }
  catch(const std::exception& e)
  {
    if (err_msg) *err_msg = DMP_pp_fun_ + e.what();
    return false;
  }
}

// ============= Load Params ================

void DMP_pp::loadParams(const std::string &filename, const std::string params_prefix)
{
  try
  {
    YAML::Node nh = YAML::LoadFile(filename);

    if (params_prefix.empty()) loadParams(nh);
    else
    {
      YAML::Node params_nh;
      if ( !YAML::getParam(nh, params_prefix, params_nh) ) throw std::runtime_error("Failed to load param \"" + params_prefix + "\"...");
      loadParams(params_nh);
    }
  }
  catch (std::exception &e)
  { throw std::runtime_error(DMP_pp_fun_ + e.what()); }
}

void DMP_pp::loadParams(const YAML::Node &n)
{
  if ( !YAML::getParam(n, "recursive_update", recursive_update) ) throw std::runtime_error(DMP_pp_fun_ + "Failed to load param \"recursive_update\"...");
  if ( !YAML::getParam(n, "sync_update", sync_update) ) throw std::runtime_error(DMP_pp_fun_ + "Failed to load param \"sync_update\"...");
  if ( !YAML::getParam(n, "adapt_to_robot_state", adapt_to_robot_state) ) throw std::runtime_error(DMP_pp_fun_ + "Failed to load param \"adapt_to_robot_state\"...");
  if ( !YAML::getParam(n, "Sigma_w_type", Sigma_w_type) ) throw std::runtime_error(DMP_pp_fun_ + "Failed to load param \"Sigma_w_type\"...");
  if ( !YAML::getParam(n, "Sigma0_tol", Sigma0_tol) ) throw std::runtime_error(DMP_pp_fun_ + "Failed to load param \"Sigma0_tol\"...");
  if ( !YAML::getParam(n, "enforce_cont_accel", enforce_cont_accel) ) throw std::runtime_error(DMP_pp_fun_ + "Failed to load param \"enforce_cont_accel\"...");
  
  if ( !YAML::getParam(n, "r0", r0) ) throw std::runtime_error(DMP_pp_fun_ + "Failed to load param \"r0\"...");
  if ( !YAML::getParam(n, "rf", rf) ) throw std::runtime_error(DMP_pp_fun_ + "Failed to load param \"rf\"...");
  if ( !YAML::getParam(n, "r1", r1) ) throw std::runtime_error(DMP_pp_fun_ + "Failed to load param \"r1\"...");
  if ( !YAML::getParam(n, "rv", rv) ) throw std::runtime_error(DMP_pp_fun_ + "Failed to load param \"rv\"...");

  if (Sigma_w_type.compare("pos") && Sigma_w_type.compare("vel") && Sigma_w_type.compare("accel") && Sigma_w_type.compare("exp"))
    throw std::runtime_error(DMP_pp_fun_ + "Unrecognized \"Sigma_w_type\":" + Sigma_w_type + "...");
}