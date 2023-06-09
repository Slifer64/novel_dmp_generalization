#ifndef $_PROJECT_384$_ONLINE_ADAPT_DMP_pp_H
#define $_PROJECT_384$_ONLINE_ADAPT_DMP_pp_H

#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <armadillo>

#include <yaml-cpp/yaml.h>

#include <gmp_lib/GMP/GMP.h>
#include <gmp_lib/GMP/GMP_Update.h>
#include <gmp_lib/GMP/GMPo.h>
#include <gmp_lib/GMP/GMPo_Update.h>

using namespace as64_;

class DMP_pp
{
public:

  DMP_pp();
  ~DMP_pp();

  void reset(bool reset_model=true);

  arma::vec getRefPos(double s) const;
  arma::vec getRefVel(double s, double s_dot) const;
  arma::vec getRefAccel(double s, double s_dot, double s_ddot) const;

  arma::vec getRefQuat(double s) const;
  arma::vec getRefRotVel(double s, double s_dot) const;
  arma::vec getRefRotAccel(double s, double s_dot, double s_ddot) const;

  void initUpdate(double s_start, const arma::vec &P0, const arma::vec &Pg, const arma::vec &Q0, const arma::vec &Qg, double Tf);

  void updatePos(const arma::vec &Pg, double s, double s_dot, const arma::vec &P, const arma::vec &P_dot);
  void updateOrient(const arma::vec &Qg, double s, double s_dot, const arma::vec &Q, const arma::vec &vRot);

  double updateViapoint(double s_v, const arma::vec &pos, const arma::vec &quat, bool reverse=false);

  std::vector<double> updateViapoints(double s, const std::vector<arma::vec> &Pv, const std::vector<arma::vec> &Qv, const std::string &vp_set, std::vector<arma::vec> *Pv2=NULL, std::vector<arma::vec> *Qv2=NULL);
  void downdateViapoints(const std::string &vp_set);

  void resetSigmaw();

  void loadParams(const std::string &filename, const std::string params_prefix="");
  void loadParams(const YAML::Node &n);

  bool save(const std::string &filename, std::string *err_msg=0);
  bool load(const std::string &filename, std::string *err_msg=0);

  arma::vec getMotionDirection(double s) const;

  bool adapt_to_robot_state;

protected:

  struct ViapointSet
  {
    ViapointSet() {}
    ViapointSet(const std::vector<double> &sv, const std::vector<arma::vec> &Pv, const std::vector<arma::vec> &Qv): s(sv), P(Pv), Q(Qv) {}
    std::vector<double> s;
    std::vector<arma::vec> P;
    std::vector<arma::vec> Q;
  };

  std::map<std::string, ViapointSet> vp_map;

  double findClosest(double s0, double sf, int n_points, const arma::vec &pos, arma::vec *p2=NULL) const;

  struct
  {
    arma::vec P0;
    arma::vec Pg;
    double s;
    double s_dot;
    bool initialized;
  } prev_pos_state;

  struct
  {
    arma::vec Q0;
    arma::vec Qg;
    double s;
    double s_dot;
    bool initialized;
  } prev_quat_state;

  // ------- Model ---------

  arma::mat W0_p;
  arma::mat W0_o;

  gmp_::GMP::Ptr gmp_p;
  gmp_::GMPo::Ptr gmp_o;

  gmp_::GMP_Update::Ptr gmp_p_up;
  gmp_::GMPo_Update::Ptr gmp_o_up;

  bool initialized;

  bool recursive_update;
  bool sync_update;
  std::string Sigma_w_type; // "pos", "vel", "accel", "exp"
  double Sigma0_tol;

  bool enforce_cont_accel;

  arma::vec r0, rf, r1;
  double rv;
};


#endif // $_PROJECT_384$_ONLINE_ADAPT_MODEL_H
