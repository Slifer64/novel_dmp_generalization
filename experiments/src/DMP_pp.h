#ifndef $_PROJECT_384$_ONLINE_ADAPT_DMP_pp_H
#define $_PROJECT_384$_ONLINE_ADAPT_DMP_pp_H

#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <armadillo>

#include <gmp_lib/GMP/GMP.h>
#include <gmp_lib/GMP/GMP_Update.h>
#include <gmp_lib/GMP/GMPo.h>
#include <gmp_lib/GMP/GMPo_Update.h>

using namespace as64_;

// DMP++
class DMP_pp
{
public:

  /** Initializes the hyper-paremeters
   */ 
  DMP_pp();

  /** Reinitializes the model for a new execution
   * @param[in] reset_model to reset the DMP weights to the initial weights learned from the demo (optional, default=true)
   */ 
  void reset(bool reset_model=true);
  
  /** Returns the DMP++ forcing term Cartesian position, velocity and acceleration.
   * @param[in] s phase variable
   * @param[in] s_dot phase variable 1st time derivative
   * @param[in] s_ddot phase variable 2nd time derivative
   */ 
  arma::vec getRefPos(double s) const;
  arma::vec getRefVel(double s, double s_dot) const;
  arma::vec getRefAccel(double s, double s_dot, double s_ddot) const;

  /** Returns the DMP++ forcing term orientation (as unit quaternion), rotational velocity and acceleration.
   * @param[in] s phase variable
   * @param[in] s_dot phase variable 1st time derivative
   * @param[in] s_ddot phase variable 2nd time derivative
   */ 
  arma::vec getRefQuat(double s) const;
  arma::vec getRefRotVel(double s, double s_dot) const;
  arma::vec getRefRotAccel(double s, double s_dot, double s_ddot) const;

  /** Initializes the weights optimization (i=0)
   * @param[in] s_start the initial value of the phase variable
   * @param[in] P0 initial Cartesian position
   * @param[in] Pg target Cartesian position
   * @param[in] Q0 initial Cartesian orientation (as unit quat)
   * @param[in] Qg target Cartesian orientation (as unit quat)
   * @param[in] Tf motion duration
   */
  void initUpdate(double s_start, const arma::vec &P0, const arma::vec &Pg, const arma::vec &Q0, const arma::vec &Qg, double Tf);

  /** Updates the position DMP weights
   * @param[in] Pg target Cartesian position
   * @param[in] s current value of phase variable
   * @param[in] s_dot current value of phase variable 1st time derivative
   * @param[in] P current Cartesian position
   * @param[in] P_dot current Cartesian velocity
   */
  void updatePos(const arma::vec &Pg, double s, double s_dot, const arma::vec &P, const arma::vec &P_dot);

  /** Updates the orientation DMP weights
   * @param[in] Qg target Cartesian orientation (as unit quat)
   * @param[in] s current value of phase variable
   * @param[in] s_dot current value of phase variable 1st time derivative
   * @param[in] Q current Cartesian orientation
   * @param[in] vRot current rotational velocity
   */
  void updateOrient(const arma::vec &Qg, double s, double s_dot, const arma::vec &Q, const arma::vec &vRot);

  /** Updates the DMP weights given a via-point
   * @param[in] s_v the initial phase variable value from which to start looking for s_min = argmin_{s}|| pos - y_s(s) ||
   * @param[in] pos via-point Cartesian position
   * @param[in] quat via-point Cartesian orientation
   * @param[in] reverse whether to look in the forward or reverse direction for s_min
   * @return s_min
   */
  double updateViapoint(double s_v, const arma::vec &pos, const arma::vec &quat, bool reverse=false);

  /** Updates the DMP weights given a set of via-points
   * @param[in] s the initial phase variable value from which to start looking for s_min = argmin_{s}|| pos - y_s(s) ||
   * @param[in] Pv std::vector with the via-point Cartesian positions
   * @param[in] quat std::vector with the via-point Cartesian orientations
   * @param[in] vp_set string used as an identifier for this set of via-points.
   * @return std::vector of the found phase variable values corresponding to each via-point
   */
  std::vector<double> updateViapoints(double s, const std::vector<arma::vec> &Pv, const std::vector<arma::vec> &Qv, const std::string &vp_set);

  /** Updates the DMP weights by removing a previous set of via-points
   * @param[in] vp_set string used as an identifier for this set of via-points.
   */
  void downdateViapoints(const std::string &vp_set);

  /** Reinitializes the "P" matrices
   */
  void resetSigmaw();

  /** Saves the position and orientation DMPs to the disk.
   * @param[in] filename the filename where to store the DMP models.
   * @param[in] err_msg Used to return an error message on failure.
   * @return true on sucess, otherwise false
   */
  bool save(const std::string &filename, std::string *err_msg=0);

  /** Loads the position and orientation DMPs from the disk.
   * @param[in] filename the filename that containts the DMP models.
   * @param[in] err_msg Used to return an error message on failure.
   * @return true on sucess, otherwise false
   */
  bool load(const std::string &filename, std::string *err_msg=0);

  bool adapt_to_robot_state; ///< flag to adapt to the current robot state

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


  /** Searches for s_min = argmin_{s}|| pos - y_s(s) || for s uniformly samples in [s0, sf].
   * @return s_min
   */
  double findClosest(double s0, double sf, int n_points, const arma::vec &pos) const;

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

  // initial DMP weights for pos and orient, learned from a demo
  arma::mat W0_p;
  arma::mat W0_o;

  // DMP for position and orientation
  gmp_::GMP::Ptr gmp_p;
  gmp_::GMPo::Ptr gmp_o;

  // Objects for updating the DMP weights
  gmp_::GMP_Update::Ptr gmp_p_up;
  gmp_::GMPo_Update::Ptr gmp_o_up;

  bool initialized;

  bool recursive_update;
  bool sync_update;
  std::string Sigma_w_type; // to choose which cost function to optimize: "pos", "vel", "accel", "exp"
  double Sigma0_tol; // for better numerical stability

  bool enforce_cont_accel; // to use also the acceleration y_ddot{i-1} in the current state constraint

  // epsilon values used in the "R" matrix
  arma::vec r0; // for initial position
  arma::vec rf; // for target position
  arma::vec r1; // for current state
  double rv; // for via-point
};


#endif // $_PROJECT_384$_ONLINE_ADAPT_MODEL_H
