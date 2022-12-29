#include <vector>
#include <DMP_pp.h>
#include <robot/robot.h>


class Controller
{
public:

// inertia, damping and stiffness for position and orientation DMP 
struct 
{
  double Mp = 5;
  double Dp = 10;
  double Kp = 40;

  double Mo = 2;
  double Do = 0.6;
  double Ko = 2;
} ctrl_params;

// via-points w.r.t. the box target, for properly inserting the carton inside
std::vector<double> place_vp_offsets = {0.15, 0.12, 0.09, 0.06, 0.03};

DMP_pp pick_model;
DMP_pp place_model;

arma::vec P0, Q0; // initial position and orientation (as unit quat)

struct
{
  arma::vec P;
  arma::vec Q;
} current_pose;

/*
Other auxiliary functions and class objects, such as:
 - robot: for communicating with the robot
 - target_reader: to get the carton's pose from an april-tag
 - vp_reader: to get the via-points pose from an april-tag
 - box_reader: to get the box's pose from an april-tag
 etc.
*/

// ===============================================================
// ===============================================================

void loadModels(const std::string &pick_model_filename, const std::string &place_model_filename)
{
  bool success;
  std::string err_msg;

  success = pick_model.load(pick_model_filename, &err_msg);
  if (!success) throw std::runtime_error(err_msg);

  success = place_model.load(place_model_filename);
  if (!success) throw std::runtime_error(err_msg);
}

// ===============================================================
// ===============================================================


Result executePick(bool reverse=false)
{
  auto &model = pick_model;

  if (reverse) PRINT_INFO_MSG("Executing reverse...\n");
  else PRINT_INFO_MSG("Executing forward...\n");

  arma::vec P_current = robot->getTaskPosition();
  arma::vec Q_current = robot->getTaskOrientation();
  if (reverse)
  {
    P_current = current_pose.P;
    Q_current = current_pose.Q;
  }

  if (reverse && arma::dot(Q_current, this->Q0) < 0) Q_current = -Q_current;

  double dt = robot->getCtrlCycle();

  arma::vec P = P_current;
  arma::vec P_dot = arma::vec().zeros(3);
  arma::vec P_ddot = arma::vec().zeros(3);

  arma::vec Q = Q_current;
  arma::vec vRot = arma::vec().zeros(3);
  arma::vec vRot_dot = arma::vec().zeros(3);

  arma::vec Fext = arma::vec().zeros(6);

  arma::vec Pg;
  arma::vec Qg;
  arma::vec Qg_prev;

  double Tf = 10; // set the pick duration to 10 sec

  double sd_dot = 1/Tf;
  double s_start = 0;
  double s_end = 1;

  if (!reverse)
  {
    this->P0 = P_current;
    this->Q0 = Q_current;
  }

  if (reverse)
  {
    Pg = P_current;
    Qg = Q_current;
    Qg_prev = Qg;

    sd_dot = -sd_dot;
    s_start = 1;
    s_end = 0;
  }

  double t = 0;
  double s = s_start;
  double s_dot = sd_dot;
  double s_ddot = 0;

  if (!reverse)
  {
    // get target for the first time
    PRINT_INFO_MSG("Waiting to receive target pose...\n");
    while (run_)
    {
      bool got_target = target_reader->getPose(&Pg, &Qg);
      if (!got_target) return Result(false, "Failed to read the target pose...");
      else break;
    }
    PRINT_SUCCESS_MSG("Got target pose!\n");

    if (arma::dot(Qg, Q) < 0) Qg = -Qg;
    Qg_prev = Qg;
  }
  
  PRINT_INFO_MSG("Executing controller...\n");

  Result result;

  double s_v = reverse ? 1.0 : 0.0;
  arma::vec P_vp = {0, 0, 0};
  arma::vec Q_vp = {1, 0, 0, 0};
  bool processed_vp = false;

  if (!reverse && enable_model_adapt()) model.resetSigmaw();

  // initialize model
  model.initUpdate(s_start, P0, Pg, Q0, Qg, Tf);

  Fext = getTaskWrench(Fext);

  while (run_)
  {
    result = updateRobot();
    if (!result.success) return result;

    Fext = getTaskWrench(Fext);
    arma::vec Fp = Fext.subvec(0, 2);
    arma::vec Fo = Fext.subvec(3, 5);

    if (!reverse) // read target only in forward execution, for retraction it is not needed
    {
      bool got_target = target_reader->getPose(&Pg, &Qg);
      if (!got_target) return Result(false, "Failed to read the target pose...");
      if (arma::dot(Qg, Qg_prev) < 0) Qg = -Qg;
      Qg_prev = Qg;
    }

    // canonical system
    double s1_dot = sd_dot;
    if (phase_stop())
    {
      double F_norm = arma::norm(Fp);
      if (arma::norm(P - Pg) < 0.03) F_norm = std::max(F_norm - 3, 0.0); // to allow some fixed contact force to be applied at the end
      else F_norm = std::max(F_norm - 0.5, 0.0); // remove small noise

      s1_dot = sd_dot / (1 + 10*arma::norm(Fp));
    }
    s_ddot = 30*(s1_dot - s_dot);

    // s_v should be not be behind the current s in forward and accordingly in reverse.
    if (reverse) s_v = std::min(s_v, s);
    else s_v = std::max(s_v, s);
    
    // allow processing of via-points both in forward and reverse
    bool got_vp = vp_reader->getPose(&P_vp, &Q_vp);

    model.adapt_to_robot_state = enable_model_adapt();
    if (!reverse && arma::norm(P - Pg) < 0.02) model.adapt_to_robot_state = false; // if too close to the target don't adapt to contact forces
    double vp_dist_thres = 0.7; // consider only via-points within this distance from the current position
    if (!processed_vp && got_vp && arma::norm(P_vp - P) < vp_dist_thres)
    {
      if (arma::dot(Q_vp, Q) < 0) Q_vp = -Q_vp;
      s_v = model.updateViapoint(s_v, P_vp, Q_vp, reverse);
    }
    model.updatePos(Pg, s, sd_dot, P, P_dot);
    model.updateOrient(Qg, s, sd_dot, Q, vRot);

    // ------ DMP++ generalized reference --------
    arma::vec Pd = model.getRefPos(s);
    arma::vec Pd_dot = model.getRefVel(s, s_dot);
    arma::vec Pd_ddot = model.getRefAccel(s, s_dot, s_ddot);
    
    // feed-forward gain term to subdue the effect of the feed-forward acceleration in the presence of external forces
    double ff = 1.0;
    ff = 1 - std::pow(std::min(arma::norm(Fp)*2, 1.0), 0.5);

    // Cartesian position DMP equation
    P_ddot = ff*Pd_ddot + ctrl_params.Dp*(Pd_dot - P_dot) + ctrl_params.Kp*(Pd - P) + Fp/ctrl_params.Mp;

    // ------ Orientation ctrl --------
    arma::vec Qd = model.getRefQuat(s);
    if (arma::dot(Qd, Q) < 0) Qd = -Qd;
    arma::vec vRotd =  model.getRefRotVel(s, s_dot);
    arma::vec vRotd_dot =  model.getRefRotAccel(s, s_dot, s_ddot);
    
    // feed-forward gain term to subdue the effect of the feed-forward acceleration in the presence of external torques
    ff = 1 - std::pow(std::min(arma::norm(Fo)*2, 1.0), 0.5);

    // orientation DMP equation
    vRot_dot = ff*vRotd_dot + ctrl_params.Do*(vRotd - vRot) + ctrl_params.Ko*math_::quatLog(math_::quatDiff(Qd, Q)) + Fo/ctrl_params.Mo;

    arma::vec V_cmd = arma::join_vert(P_dot, vRot);
    robot->setTaskVelocity(V_cmd, P, Q);
    
    // numerical integration
    t = t + dt;
    
    s = s + s_dot*dt;
    s_dot = s_dot + s_ddot*dt;

    P = P + P_dot*dt;
    P_dot = P_dot + P_ddot*dt;

    Q = math_::quatProd(math_::quatExp(vRot*dt), Q);
    vRot = vRot + vRot_dot*dt;

    bool s_finished = reverse ? s < s_end : s > s_end;

    if (s_finished) break;
  }

  if (!reverse) gripper->close(); // close gripper only in forward

  robot->setTaskVelocity(arma::vec().zeros(6));

  PRINT_SUCCESS_MSG("Finished execution!\n");

  return Result(true);
}


// ===============================================================
// ===============================================================


Result executePlace(bool reverse=false)
{
  auto &model = place_model;

  // enable phase stopping. This also automatically disables adaptation to the current robot state
  gui->triggerPhaseStop();

  if (reverse) PRINT_INFO_MSG("Retract from place...\n");
  else PRINT_INFO_MSG("Executing place...\n");

  arma::vec P_current = robot->getTaskPosition();
  arma::vec Q_current = robot->getTaskOrientation();

  if (reverse && arma::dot(Q_current, this->Q0) < 0) Q_current = -Q_current;

  double dt = robot->getCtrlCycle();

  arma::vec P = P_current;
  arma::vec P_dot = arma::vec().zeros(3);
  arma::vec P_ddot = arma::vec().zeros(3);

  arma::vec Q = Q_current;
  arma::vec vRot = arma::vec().zeros(3);
  arma::vec vRot_dot = arma::vec().zeros(3);

  arma::vec Fext = arma::vec().zeros(6);

  arma::vec Pg;
  arma::vec Qg;
  arma::vec Qg_prev;

  double Tf = 8;  // motion duration for place
  double sd_dot = 1/Tf;
  double s_start = 0;
  double s_end = 1;

  if (!reverse)
  {
    this->P0 = P_current;
    this->Q0 = Q_current;
  }

  if (reverse)
  {
    Pg = P_current;
    Qg = Q_current;
    Qg_prev = Qg;

    sd_dot = -sd_dot;
    s_start = 1;
    s_end = 0;
  }

  double t = 0;

  double s = s_start;
  double s_dot = sd_dot;
  double s_ddot = 0;

  if (!reverse)
  {
    // get target for the first time
    PRINT_INFO_MSG("Waiting to receive target pose...\n");
    while (run_)
    {
      int got_target = box_reader->getPose(&Pg, &Qg);
      if (got_target < 0) return Result(false, "Failed to read the box pose...");
      if (got_target == 0) break;
    }
    PRINT_SUCCESS_MSG("Got box pose!\n");

    if (arma::dot(Qg, Q) < 0) Qg = -Qg;
    Qg_prev = Qg;
  }
  
  PRINT_INFO_MSG("Executing controller...\n");

  Result result;

  // initialize model
  std::vector<arma::vec> vp_offsets;
  for (double offset : place_vp_offsets) vp_offsets.push_back({0, 0, -offset});
  std::vector<arma::vec> Pv, Qv;
  if (!reverse)
  {
    model.reset();
    model.initUpdate(s_start, P0, Pg, Q0, Qg, Tf);
    model.adapt_to_robot_state = false;

    // express the via-points in the robot's base frame
    arma::mat Rg = math_::quat2rotm(Qg);
    for (int i=0; i<vp_offsets.size(); i++)
    {
      Pv.push_back(Pg + Rg*vp_offsets[i]);
      Qv.push_back(Qg);
    }

    model.updateViapoints(0, Pv, Qv, "place-vp");
  }
  else
  {
    model.initUpdate(s_start, P0, Pg, Q0, Qg, Tf);
  }

  while (run_)
  {
    result = updateRobot();
    if (!result.success) return result;

    Fext = getTaskWrench(Fext);
    arma::vec Fp = Fext.subvec(0, 2);
    arma::vec Fo = Fext.subvec(3, 5);

    double s1_dot = sd_dot;
    if (phase_stop()) s1_dot = sd_dot / (1 + 10*arma::norm(Fp));
    s_ddot = 30*(s1_dot - s_dot);

    if (!reverse && box_reader->getPose(&Pg, &Qg)) // got_new_target
    {
      Pv.clear();
      Qv.clear();
  
      // phase stopping will change s_dot
      // we don't want to adapt the trajectory to the current state, so use sd_dot instead
      model.updatePos(Pg, s, sd_dot, P, P_dot); 
      model.updateOrient(Qg, s, sd_dot, Q, vRot);

      model.downdateViapoints("place-vp");
      // calc new place via-points
      arma::mat Rg = math_::quat2rotm(Qg);
      Pv.reserve(vp_offsets.size());
      Qv.reserve(vp_offsets.size());
      for (int i=0; i<vp_offsets.size(); i++)
      {
        Pv.push_back(Pg + Rg*vp_offsets[i]);
        Qv.push_back(Qg);
      }
      // model.updateViapoints(s, Pv, Qv, "place-vp");
      model.updateViapoints(0, Pv, Qv, "place-vp");
    }

    // ------ DMP++ generalized reference --------
    arma::vec Pd = model.getRefPos(s);
    arma::vec Pd_dot = model.getRefVel(s, s_dot);
    arma::vec Pd_ddot = model.getRefAccel(s, s_dot, s_ddot);
    
    // position DMP equation
    P_ddot = Pd_ddot + ctrl_params.Dp*(Pd_dot - P_dot) + ctrl_params.Kp*(Pd - P) + Fp/ctrl_params.Mp;

    // ------ Orientation ctrl --------
    arma::vec Qd = model.getRefQuat(s);
    if (arma::dot(Qd, Q) < 0) Qd = -Qd;
    arma::vec vRotd =  model.getRefRotVel(s, s_dot);
    arma::vec vRotd_dot =  model.getRefRotAccel(s, s_dot, s_ddot);
    
    // orientation DMP equation
    vRot_dot = vRotd_dot + ctrl_params.Do*(vRotd - vRot) + ctrl_params.Ko*math_::quatLog(math_::quatDiff(Qd, Q)) + Fo/ctrl_params.Mo;

    arma::vec V_cmd = arma::join_vert(P_dot, vRot);
    robot->setTaskVelocity(V_cmd, P, Q);
    
    // numerical integration
    t = t + dt;
    
    s = s + s_dot*dt;
    s_dot = s_dot + s_ddot*dt;

    P = P + P_dot*dt;
    P_dot = P_dot + P_ddot*dt;

    Q = math_::quatProd(math_::quatExp(vRot*dt), Q);
    vRot = vRot + vRot_dot*dt;

    bool s_finished = reverse ? s < s_end : s > s_end;
    if (s_finished) break;
  }

  if (!reverse) gripper->open();

  robot->setTaskVelocity(arma::vec().zeros(6));

  PRINT_SUCCESS_MSG("Finished execution!\n");

  return Result(true);
}


} // Controller