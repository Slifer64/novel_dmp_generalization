#include <online_adapt_controller/controller.h>

#include <iostream>
// #include <io_lib/xml_parser.h>
#include <io_lib/file_io.h>
#include <plot_lib/qt_plot.h>
#include <math_lib/math_lib.h>

#include <ros/ros.h>
#include <ros/package.h>
#include <tf/transform_broadcaster.h>

#include <yaml-cpp/yaml.h>

#include <gmp_lib/io/gmp_io.h>

using namespace as64_;

#define OnlineAdaptController_fun_ std::string("[OnlineAdaptController::") + __func__ + "]: "


std::vector<std::string> LogState::var_names = {"t", "s", "P", "V", "V_dot", "Pd", "Vd", "Vd_dot", "Pg", "Fext", "sv", "Pv", "obst"};

// ================================================================
// ================================================================

OnlineAdaptController::OnlineAdaptController(MainController *main_ctrl, const std::string &ctrl_name)
: Controller(main_ctrl, ctrl_name)
{
  this->main_ctrl = main_ctrl;
  this->robot = main_ctrl->robot.get();

  ros::NodeHandle nh("~");

  if (!nh.getParam("online_adapt_controller_data_path", default_data_path)) default_data_path = "/config/";
  default_data_path = ros::package::getPath("online_adapt_controller") + "/" + default_data_path;

  std::string robot_descr_param;
  if (!nh.getParam("robot_description_name", robot_descr_param)) throw std::runtime_error("\33[1m\33[31mFailed to load param 'robot_description_name'...\33[0m\n");
  if (!nh.getParam("base_link", base_link)) throw std::runtime_error("\33[1m\33[31mFailed to load param 'base_link'...\33[0m\n");
  std::string tool_link;
  if (!nh.getParam("tool_link", tool_link)) throw std::runtime_error("\33[1m\33[31mFailed to load param 'tool_link'...\33[0m\n");

  std::cerr << "==========> Ok 132\n";

  std::string err_msg;
  if (!reloadAprilTagListener(&err_msg)) throw std::runtime_error(err_msg);

  this->marker_array_topic = "/" + getNodeNameID() + "_online_adapt_ctrl_marker_topic";

  viz.reset(new online_adapt_::Visualizer(marker_array_topic, base_link));

  gripper.reset(new online_adapt_::Gripper("/epick"));

  pick_target_reader.reset(new online_adapt_::TagReader());

  getPickTargetPose = [this](arma::vec *P, arma::vec *Q){ return pick_target_reader->getPose(P, Q); };

  vp_reader.reset(new online_adapt_::TagReader());
  place_target_reader.reset(new online_adapt_::TagReader());
  obstacle_reader.reset(new online_adapt_::TagReader());

  std::cerr << "==========> Ok 168\n";

  ExecResultMsg msg = loadParams();
  if (msg.getType() != ExecResultMsg::INFO) online_adapt_::PRINT_ERROR_MSG(msg.getMsg());

  bool success;
  success = pick_model.load(pick_model_filename, &err_msg);
  if (!success) throw std::runtime_error(err_msg);
  success = place_model.load(place_model_filename);
  if (!success) throw std::runtime_error(err_msg);

  std::cerr << "==========> Ok 265\n";

  pick_target_reader->setTagReadFunction([this](arma::vec *p, arma::vec *Q){ return this->getTagPose(pick_target_reader->tagId(), p, Q); });
  vp_reader->setTagReadFunction([this](arma::vec *p, arma::vec *Q){ return this->getTagPose(vp_reader->tagId(), p, Q); });
  place_target_reader->setTagReadFunction([this](arma::vec *p, arma::vec *Q){ return this->getTagPose(place_target_reader->tagId(), p, Q); });
  obstacle_reader->setTagReadFunction([this](arma::vec *p, arma::vec *Q){ return this->getTagPose(obstacle_reader->tagId(), p, Q); });

  std::cerr << "==========> Ok 347\n";

  // launchTagCtrlSignalsThread();

  launchExternalGripperCtrlThread();
}

OnlineAdaptController::~OnlineAdaptController()
{
  keep_alive = false;
  if (tag_ctrl_thead.joinable()) tag_ctrl_thead.join();
}

void OnlineAdaptController::launchExternalGripperCtrlThread()
{
  std::thread([this]()
  {
    while (keep_alive)
    {
      if (this->main_ctrl->getCustomDigitalIn())
      {
        if (gripper->isOpen()) gripper->close();
        else gripper->open();
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }
  }).detach();
}

// ===============================================================

void OnlineAdaptController::launchTagCtrlSignalsThread()
{
  tag_ctrl_thead = std::thread([this]()
  {
    while (keep_alive)
    {
      if (tag_listener->getTagDetection(stop_exec_tag_id, true).tag_id > 0) gui->triggerStopExec();
      else if (tag_listener->getTagDetection(phase_stop_tag_id, true).tag_id > 0) gui->triggerPhaseStop();
      else if (tag_listener->getTagDetection(adapt_to_robot_tag_id, true).tag_id > 0) gui->triggerAdaptToRobot();
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
  });
}

bool OnlineAdaptController::getTagPose(int tag_id, arma::vec *P, arma::vec *Q) const
{
  apriltag_ros::AprilTagListener::TagDetection tag = tag_listener->getTagDetection(tag_id);

  if (tag.tag_id < 0) return false;

  *P = { tag.pos.x(), tag.pos.y(), tag.pos.z() };
  *Q = {tag.quat.w(), tag.quat.x(), tag.quat.y(), tag.quat.z()};

  return true;
}

QPushButton *OnlineAdaptController::createGui(MainWindow *parent)
{
  QPushButton *btn = new QPushButton(this->ctrl_name.c_str());
  this->gui = new OnlineAdaptWin(this, parent);
  QObject::connect( btn, &QPushButton::pressed, parent, [this](){ gui->launch(); });
  return btn;
}

// ===============================================================

void OnlineAdaptController::start()
{
  Result result;
  ExecResultMsg msg;
  
  // ------ Load params -------
  msg = loadParams();
  if (msg.getType() != ExecResultMsg::INFO)
  {
    emit gui->stopCtrlSignal(msg);
    exec_stop_sem.notify();
    return;
  }

  this->setMode(rw_::CART_VEL_CTRL);
  run_ = true;

  // ------ Init execution -------

  result = init();
  if (!result.success)
  {
    emit gui->stopCtrlSignal(ExecResultMsg(ExecResultMsg::ERROR, result.err_msg));
    exec_stop_sem.notify();
    return;
  }

  if (!skip_pick)
  {
    // ------ Forward execution -------
    std::thread thr = std::thread([this, &result](){ result = executePick(false); });
    int err_code = thr_::setThreadPriority(thr, SCHED_FIFO, 99);
    if (err_code) online_adapt_::PRINT_WARNING_MSG("[OnlineAdaptController::start]: Failed to set thread priority! Reason:\n" + thr_::setThreadPriorErrMsg(err_code));
    if (thr.joinable()) thr.join();
    if (!result.success)
    {
      emit gui->stopCtrlSignal(ExecResultMsg(ExecResultMsg::ERROR, result.err_msg));
      exec_stop_sem.notify();
      return;
    }

    // Wait
    // while (run_ && !execute_reverse) std::this_thread::sleep_for(std::chrono::milliseconds(30));

    // ------ Reverse execution -------
    if (run_ && execute_reverse)
    {
      Result result = executePick(true);
      if (!result.success)
      {
        emit gui->stopCtrlSignal(ExecResultMsg(ExecResultMsg::ERROR, result.err_msg));
        exec_stop_sem.notify();
        return;
      }
    }
  }
  
  // Wait
  while (run_ && !execute_place) std::this_thread::sleep_for(std::chrono::milliseconds(30));

  // ------ Place execution -------
  if (run_ && execute_place)
  {
    Result result = executePlace(false);
    if (!result.success)
    {
      emit gui->stopCtrlSignal(ExecResultMsg(ExecResultMsg::ERROR, result.err_msg));
      exec_stop_sem.notify();
      return;
    }
  }

  // Wait
  while (run_ && !retract_from_place) std::this_thread::sleep_for(std::chrono::milliseconds(30));

  // ------ Retract from place -------
  if (run_ && retract_from_place)
  {
    Result result = executePlace(true);
    if (!result.success)
    {
      emit gui->stopCtrlSignal(ExecResultMsg(ExecResultMsg::ERROR, result.err_msg));
      exec_stop_sem.notify();
      return;
    }
  }

  // ------ Wait for termination -------
  online_adapt_::PRINT_INFO_MSG("Finished. Waiting command from GUI...\n");
  while (run_) std::this_thread::sleep_for(std::chrono::milliseconds(30));

  online_adapt_::PRINT_SUCCESS_MSG("Got termination signal!\n");

  exec_stop_sem.notify();
}

ExecResultMsg OnlineAdaptController::stop()
{
  run_ = false;
  ExecResultMsg msg;

  if ( exec_stop_sem.wait_for(1500) ) msg = ExecResultMsg(ExecResultMsg::INFO, "Controller stopped!");
  else msg = ExecResultMsg(ExecResultMsg::WARNING, "Time-out reached on waiting for controller to stop...");

  robot->setTaskVelocity(arma::vec().zeros(6));

  setMode(rw_::IDLE);

  gripper->open();

  if (auto_save)
  {
    ExecResultMsg msg = saveLoggedData();
    if (msg.getType() != ExecResultMsg::INFO) emit gui->showMsgSignal(msg);
  }

  return msg;
}

Result OnlineAdaptController::init()
{
  try
  {
    sim = gui->runInSimulation();

    pick_target_reader->setTagReadFunction([this](arma::vec *p, arma::vec *Q){ return this->getTagPose(pick_target_reader->tagId(), p, Q); });
    vp_reader->setTagReadFunction([this](arma::vec *p, arma::vec *Q){ return this->getTagPose(vp_reader->tagId(), p, Q); });
    place_target_reader->setTagReadFunction([this](arma::vec *p, arma::vec *Q){ return this->getTagPose(place_target_reader->tagId(), p, Q); });
    obstacle_reader->setTagReadFunction([this](arma::vec *p, arma::vec *Q){ return this->getTagPose(obstacle_reader->tagId(), p, Q); });

    pick_model.reset(reset_model());
    place_model.reset(true);

    viz->stop();
    viz->start();

    x_progress = 0;

    gripper->open();

    frw_logger.clear();
    frw_logger.reserve(log_reserve);

    rev_logger.clear();
    rev_logger.reserve(log_reserve);

    frw_place_logger.clear();
    frw_place_logger.reserve(log_reserve);

    rev_place_logger.clear();
    rev_place_logger.reserve(log_reserve);

    return Result(true, "");
  }
  catch (std::exception &e)
  { return Result(false, OnlineAdaptController_fun_ + e.what()); }
}

Result OnlineAdaptController::executePick(bool reverse)
{
  auto &model = pick_model;

  if (reverse) online_adapt_::PRINT_INFO_MSG("Executing reverse...\n");
  else online_adapt_::PRINT_INFO_MSG("Executing forward...\n");

  arma::vec P_current = robot->getTaskPosition();
  arma::vec Q_current = robot->getTaskOrientation();
  if (sim && reverse)
  {
    P_current = current_pose.P;
    Q_current = current_pose.Q;
  }

  if (reverse && arma::dot(Q_current, this->Q0) < 0) Q_current = -Q_current;

  auto &logger = reverse ? rev_logger : frw_logger; 

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

  double Tf = this->pick_Tf;

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

    viz->publishFrames(false);
    viz->clearViapoints();
  }

  double t = 0;
  int count = 0;

  double s = s_start;
  double s_dot = sd_dot;
  double s_ddot = 0;

  current_pose.P = P;
  current_pose.Q = Q;

  if (!reverse)
  {
    // get target for the first time
    online_adapt_::PRINT_INFO_MSG("Waiting to receive target pose...\n");
    while (run_)
    {
      int got_target = getPickTargetPose(&Pg, &Qg);
      if (got_target == 0) return Result(false, "Failed to read the target pose...");
      if (got_target > 0) break;
    }
    online_adapt_::PRINT_SUCCESS_MSG("Got target pose!\n");

    if (arma::dot(Qg, Q) < 0) Qg = -Qg;
    Qg_prev = Qg;
  }
  
  online_adapt_::PRINT_INFO_MSG("Executing controller...\n");

  Result result;

  arma::wall_clock timer;
  std::vector<double> elaps_times;
  elaps_times.reserve(10000);

  double s_v = reverse ? 1.0 : 0.0;
  arma::vec P_vp = {0, 0, 0};
  arma::vec Q_vp = {1, 0, 0, 0};
  bool processed_vp = false;

  int log_count = log_every;

  if (!reverse && enable_model_adapt()) model.resetSigmaw();

  // initialize model
  model.initUpdate(s_start, P0, Pg, Q0, Qg, Tf);

  Fext = getTaskWrench(Fext);
  if (!reverse && arma::norm(Fext) > 1e-2)
  {
    return Result(false, "Non-zero wrench detected!\n||Fext|| = " + std::to_string(arma::norm(Fext)) +"\nMake sure the F/T sensor is properly biased.");
  }

  double d_phase_stop_prev = 0;
  double fv_prev = 0;

  arma::vec obst_pose = arma::vec().ones(7) * std::nan("");;

  while (run_)
  {
    if (reverse) x_progress = std::max(s, 0.0);
    else x_progress = std::min(s, 1.0);

    result = updateRobot();
    if (!result.success) return result;

    Fext = getTaskWrench(Fext);
    arma::vec Fp = Fext.subvec(0, 2);
    arma::vec Fo = Fext.subvec(3, 5);

    if (!reverse) // read target only in forward execution
    {
      int got_target = getPickTargetPose(&Pg, &Qg);
      if (got_target < 0) return Result(false, "Failed to read the target pose...");
      if (arma::dot(Qg, Qg_prev) < 0) Qg = -Qg;
      Qg_prev = Qg;
    }

    double s1_dot = sd_dot;
    double fv = 0;
    if (phase_stop())
    {
      double F_norm = arma::norm(Fp);
      if (arma::norm(P - Pg) < 0.03) F_norm = std::max(F_norm - 3, 0.0); // to allow some fixed contact force to be applied at the end
      else F_norm = std::max(F_norm - 0.5, 0.0); // remove small noise

      double a1 = 0.05;
      double s_stop = 10;
      double d = (1-a1)*d_phase_stop_prev + a1*arma::norm(Fp);
      d_phase_stop_prev = d;
      s1_dot = sd_dot / (1 + s_stop*d);
      fv = 0.1*arma::dot(Fp, model.getMotionDirection(s));
      double a_fv = 0.05;
      fv = (1-a_fv)*fv_prev + a_fv*fv;
      fv_prev = fv;
    }
    s_ddot = 30*(s1_dot - s_dot) + fv;

    // s_v should be not be behind the current s in forward and accordingly in reverse.
    if (reverse) s_v = std::min(s_v, s);
    else s_v = std::max(s_v, s);
    
    // allow processing of via-points both in forward and reverse
    bool got_vp = vp_reader->getPoseHelper(&P_vp, &Q_vp);

    timer.tic();

    model.adapt_to_robot_state = enable_model_adapt();
    if (!reverse && arma::norm(P - Pg) < 0.02) model.adapt_to_robot_state = false; // if too close to the target don't adapt to contact forces
    if (!processed_vp && got_vp && arma::norm(P_vp - P) < vp_dist_thres)
    {
      if (arma::dot(Q_vp, Q) < 0) Q_vp = -Q_vp;
      s_v = model.updateViapoint(s_v, P_vp, Q_vp, reverse);
      // processed_vp = true;
    }
    model.updatePos(Pg, s, sd_dot, P, P_dot);
    model.updateOrient(Qg, s, sd_dot, Q, vRot);

    elaps_times.push_back(timer.toc());

    // ------ Position ctrl --------
    arma::vec Pd = model.getRefPos(s);
    arma::vec Pd_dot = model.getRefVel(s, s_dot);
    arma::vec Pd_ddot = model.getRefAccel(s, s_dot, s_ddot);
    
    double ff = 1.0;
    // if (model.adapt_to_robot_state) ff = 0;
    ff = 1 - std::pow(std::min(arma::norm(Fp)*2, 1.0), 0.5);

    P_ddot = ff*Pd_ddot + ctrl_params.Dp*(Pd_dot - P_dot) + ctrl_params.Kp*(Pd - P) + Fp/ctrl_params.Mp;

    // ------ Orientation ctrl --------
    arma::vec Qd = model.getRefQuat(s);
    if (arma::dot(Qd, Q) < 0) Qd = -Qd;
    arma::vec vRotd =  model.getRefRotVel(s, s_dot);
    arma::vec vRotd_dot =  model.getRefRotAccel(s, s_dot, s_ddot);
    
    ff = 1 - std::pow(std::min(arma::norm(Fo)*2, 1.0), 0.5);
    vRot_dot = ff*vRotd_dot + ctrl_params.Do*(vRotd - vRot) + ctrl_params.Ko*math_::quatLog(math_::quatDiff(Qd, Q)) + Fo/ctrl_params.Mo;

    arma::vec V_cmd = arma::join_vert(P_dot, vRot);
    result = checkVelLimits(V_cmd);
    if (!result.success) return result;

    if (!sim) robot->setTaskVelocity(V_cmd, P, Q);

    // log data
    if (log_count == log_every)
    {
      log("t") = {t};
      log("s") = {s, s_dot, s_ddot};
      log("P") = arma::join_vert(P, Q);
      log("V") = arma::join_vert(P_dot, vRot);
      log("V_dot") = arma::join_vert(P_ddot, vRot_dot);
      log("Pd") = arma::join_vert(Pd, Qd);
      log("Vd") = arma::join_vert(Pd_dot, vRotd);
      log("Vd_dot") = arma::join_vert(Pd_ddot, vRotd_dot);
      log("Pg") = arma::join_vert(Pg, Qg);
      log("Fext") = Fext;
      if (got_vp) log("sv") = {s_v};
      else log("sv") = {-1};
      log("Pv") = arma::join_vert(P_vp, Q_vp);
      log("obst") = obst_pose;
      logger.push_back(log);
      log_count = 0;
    }
    log_count++;
    
    // numerical integration
    t = t + dt;
    
    s = s + s_dot*dt;
    s_dot = s_dot + s_ddot*dt;

    P = P + P_dot*dt;
    P_dot = P_dot + P_ddot*dt;

    Q = math_::quatProd(math_::quatExp(vRot*dt), Q);
    vRot = vRot + vRot_dot*dt;

    current_pose.P = P;
    current_pose.Q = Q;

    double ep, eo;
    
    if (reverse)
    {
      ep = arma::norm(P - P0);
      eo = arma::norm(math_::quatLog(math_::quatDiff(Q, Q0)));
    }
    else
    {
      ep = arma::norm(P - Pg);
      eo = arma::norm(math_::quatLog(math_::quatDiff(Q, Qg)));
    }
    
    viz->setRobotState(P, Q);
    viz->setRefState(Pd, Qd);
    if (got_vp) viz->setViapoint(P_vp, Q_vp);

    bool s_finished = reverse ? s < s_end : s > s_end;

    if (s_finished && ep < 0.01 && eo < 0.02)
    {
      std::cerr << "ep = " << ep << "  ,  eo = " << eo << "\n";
      break;
    }

    if (s_finished)
    {
      std::cerr << "ep = " << ep << "  ,  eo = " << eo << "\n";
    }
  }

  if (!reverse) gripper->close(); // close gripper only in forward

  robot->setTaskVelocity(arma::vec().zeros(6));

  online_adapt_::PRINT_SUCCESS_MSG("Finished execution!\n");

  arma::rowvec elaps_t_data = arma::rowvec(elaps_times) * 1000; // convert to ms
  if (!elaps_t_data.empty())
  {
    double mean_elaps_t = arma::mean(elaps_t_data);
    double std_elaps_t = arma::stddev(elaps_t_data);
    double max_elaps_t = arma::max(elaps_t_data);
    double min_elaps_t = arma::min(elaps_t_data);

    std::cerr << "======= Elapsed time (ms) ======\n";
    std::cerr << "std_range: [" << std::max(0.,mean_elaps_t-std_elaps_t) << " -  " << mean_elaps_t + std_elaps_t <<"]\n";
    std::cerr << "mean     : " << mean_elaps_t << " +/- " << std_elaps_t <<"\n";
    std::cerr << "min      : " << min_elaps_t << "\n";
    std::cerr << "max      : " << max_elaps_t << "\n";
    std::cerr << "total_time: " << t << " sec\n";
    std::cerr << "==============================\n";
  }

  return Result(true);
}

Result OnlineAdaptController::executePlace(bool reverse)
{
  auto &model = place_model;

  // gui->triggerPhaseStop();

  if (reverse) online_adapt_::PRINT_INFO_MSG("Retract from place...\n");
  else online_adapt_::PRINT_INFO_MSG("Executing place...\n");

  arma::vec P_current = robot->getTaskPosition();
  arma::vec Q_current = robot->getTaskOrientation();
  if (sim)
  {
    P_current = current_pose.P;
    Q_current = current_pose.Q;
  }

  if (reverse && arma::dot(Q_current, this->Q0) < 0) Q_current = -Q_current;

  auto &logger = reverse ? rev_place_logger : frw_place_logger; 

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

  double Tf = this->place_Tf;

  double sd_dot = 1/Tf;
  double s_start = 0;
  double s_end = 1;

  if (!reverse)
  {
    this->P0 = P_current;
    this->Q0 = Q_current;

    viz->publishFrames(true);
  }

  if (reverse)
  {
    Pg = P_current;
    Qg = Q_current;
    Qg_prev = Qg;

    sd_dot = -sd_dot;
    s_start = 1;
    s_end = 0;

    viz->publishFrames(false);
  }

  double t = 0;
  int count = 0;

  double s = s_start;
  double s_dot = sd_dot;
  double s_ddot = 0;

  current_pose.P = P;
  current_pose.Q = Q;

  if (!reverse)
  {
    // get target for the first time
    online_adapt_::PRINT_INFO_MSG("Waiting to receive target pose...\n");
    while (run_)
    {
      int got_target = place_target_reader->getPose(&Pg, &Qg);
      if (got_target > 0) break;
      if (got_target == 0) return Result(false, "Failed to read the box pose...");
    }
    online_adapt_::PRINT_SUCCESS_MSG("Got box pose!\n");

    if (arma::dot(Qg, Q) < 0) Qg = -Qg;
    Qg_prev = Qg;
  }
  
  online_adapt_::PRINT_INFO_MSG("Executing controller...\n");

  Result result;

  int log_count = log_every;

  // initialize model
  online_adapt_::ViaPoints place_vp;
  if (!reverse)
  {
    model.reset();
    model.initUpdate(s_start, P0, Pg, Q0, Qg, Tf);
    model.adapt_to_robot_state = false;
    place_vp = place_target_reader->getViaPoints(Pg, Qg);
    model.updateViapoints(0, place_vp.pos, place_vp.quat, "place-vp");

    // double sv = 0;
    // for (int i=0; i<Pv.size(); i++) sv = model.updateViapoint(sv, Pv[i], Qv[i], false);

    viz->setPlaceViapoints(place_vp.pos, place_vp.quat);
  }
  else
  {
    model.initUpdate(s_start, P0, Pg, Q0, Qg, Tf);
  }

  arma::wall_clock timer;
  std::vector<double> elaps_times;
  elaps_times.reserve(10000);

  double d_phase_stop_prev = 0;

  online_adapt_::ViaPoints viz_obst_vp;

  arma::vec obst_pose = arma::vec().ones(7) * std::nan("");;

  while (run_)
  {
    if (reverse) x_progress = std::max(s, 0.0);
    else x_progress = std::min(s, 1.0);

    result = updateRobot();
    if (!result.success) return result;

    Fext = getTaskWrench(Fext);
    arma::vec Fp = Fext.subvec(0, 2);
    arma::vec Fo = Fext.subvec(3, 5);

    double s1_dot = sd_dot;
    if (phase_stop())
    {
      double F_norm = arma::norm(Fp);
      F_norm = std::max(F_norm - 2.5, 0.0); // remove small noise

      double a1 = 0.05;
      double s_stop = 10;
      double d = (1-a1)*d_phase_stop_prev + a1*F_norm;
      d_phase_stop_prev = d;
      s1_dot = sd_dot / (1 + s_stop*d);
    }
    s_ddot = 30*(s1_dot - s_dot);

    if (s > 1)
    {
      s_ddot = 0;
      s_dot = 0;
      s = 1;
    }

    timer.tic();

    if (!reverse && place_target_reader->getPoseHelper(&Pg, &Qg)) // got_new_target
    {
      place_vp.pos.clear();
      place_vp.quat.clear();

      if (arma::dot(Qg, Qg_prev) < 0) Qg = -Qg;
      Qg_prev = Qg;
  
      // phase stopping will change s_dot
      // we don't want to adapt the trajectory to the current state, so use sd_dot instead
      model.updatePos(Pg, s, sd_dot, P, P_dot); 
      model.updateOrient(Qg, s, sd_dot, Q, vRot);

      // model.downdateViapoints("place-vp");
      // model.updateViapoints(s, Pv, Qv, "place-vp");
      place_vp = place_target_reader->getViaPoints(Pg, Qg);
      model.updateViapoints(0, place_vp.pos, place_vp.quat, "place-vp");
    }

    arma::vec P_obst;
    arma::vec Q_obst;
    online_adapt_::ViaPoints obst_vp;
    if (!reverse && obstacle_reader->getPose(&P_obst, &Q_obst))
    {
      obst_vp = obstacle_reader->getViaPoints(P_obst, Q_obst);

      obst_pose = arma::join_vert(P_obst, Q_obst);

      for (int i=0; i<obst_vp.pos.size(); i++)
      {
        obst_vp.pos[i](2) = std::nan(""); // put 'nan' at z
        obst_vp.quat[i](0) = std::nan(""); // one nan element suffices to invalidate the entrire quaternion
      }

      // run forward the model to see if it penetrates the obstacle
      bool penetrare_obst = false;

      arma::rowvec s1 = arma::linspace<arma::rowvec>(s, std::min(s+0.25, 1.0), 30);
      arma::mat R_obst = math_::quat2rotm(Q_obst);
      for (int j=0; j<s1.size(); j++)
      {
        // express the position in the obstacle's frame
        arma::vec P = R_obst.t() * (model.getRefPos(s1(j)) - P_obst);
        penetrare_obst = !arma::any(arma::abs(P) > obst_bounds);
      }

      viz->showObstacle(P_obst, Q_obst, 2*obst_bounds(0), 2*obst_bounds(1), 2*obst_bounds(2));

      if (penetrare_obst)
      {
        // model.downdateViapoints("obst-vp");
        model.updateViapoints(0, obst_vp.pos, obst_vp.quat, "obst-vp", &obst_vp.pos, &obst_vp.quat);
        viz_obst_vp = obst_vp;
      }
    }


    elaps_times.push_back(timer.toc());

    // ------ Position ctrl --------
    arma::vec Pd = model.getRefPos(s);
    arma::vec Pd_dot = model.getRefVel(s, s_dot);
    arma::vec Pd_ddot = model.getRefAccel(s, s_dot, s_ddot);
    
    P_ddot = Pd_ddot + ctrl_params.Dp*(Pd_dot - P_dot) + ctrl_params.Kp*(Pd - P) + Fp/ctrl_params.Mp;

    // ------ Orientation ctrl --------
    arma::vec Qd = model.getRefQuat(s);
    if (arma::dot(Qd, Q) < 0) Qd = -Qd;
    arma::vec vRotd =  model.getRefRotVel(s, s_dot);
    arma::vec vRotd_dot =  model.getRefRotAccel(s, s_dot, s_ddot);
    
    vRot_dot = vRotd_dot + ctrl_params.Do*(vRotd - vRot) + ctrl_params.Ko*math_::quatLog(math_::quatDiff(Qd, Q)) + Fo/ctrl_params.Mo;

    arma::vec V_cmd = arma::join_vert(P_dot, vRot);
    result = checkVelLimits(V_cmd);
    if (!result.success) return result;

    if (!sim) robot->setTaskVelocity(V_cmd, P, Q);

    // log data
    if (log_count == log_every)
    {
      log("t") = {t};
      log("s") = {s, s_dot, s_ddot};
      log("P") = arma::join_vert(P, Q);
      log("V") = arma::join_vert(P_dot, vRot);
      log("V_dot") = arma::join_vert(P_ddot, vRot_dot);
      log("Pd") = arma::join_vert(Pd, Qd);
      log("Vd") = arma::join_vert(Pd_dot, vRotd);
      log("Vd_dot") = arma::join_vert(Pd_ddot, vRotd_dot);
      log("Pg") = arma::join_vert(Pg, Qg);
      log("Fext") = Fext;
      log("obst") = obst_pose;
      logger.push_back(log);
      log_count = 0;
    }
    log_count++;
    
    // numerical integration
    t = t + dt;
    
    s = s + s_dot*dt;
    s_dot = s_dot + s_ddot*dt;

    P = P + P_dot*dt;
    P_dot = P_dot + P_ddot*dt;

    Q = math_::quatProd(math_::quatExp(vRot*dt), Q);
    vRot = vRot + vRot_dot*dt;

    current_pose.P = P;
    current_pose.Q = Q;

    double ep, eo;
    
    if (reverse)
    {
      ep = arma::norm(P - P0);
      eo = arma::norm(math_::quatLog(math_::quatDiff(Q, Q0)));
    }
    else
    {
      ep = arma::norm(P - Pg);
      eo = arma::norm(math_::quatLog(math_::quatDiff(Q, Qg)));
    }

    viz->setRobotState(P, Q);
    viz->setRefState(Pd, Qd);
    viz->setPlaceViapoints(place_vp.pos, place_vp.quat);
    viz->setObstViapoints(viz_obst_vp.pos, viz_obst_vp.quat);

    bool s_finished = reverse ? s < s_end : s > s_end;

    if (s_finished && ep < 0.01 && eo < 0.02)
    {
      std::cerr << "ep = " << ep << "  ,  eo = " << eo << "\n";
      break;
    }

    if (s_finished)
    {
      std::cerr << "ep = " << ep << "  ,  eo = " << eo << "\n";
    }
  }

  if (!reverse) gripper->open();

  robot->setTaskVelocity(arma::vec().zeros(6));

  online_adapt_::PRINT_SUCCESS_MSG("Finished execution!\n");

  arma::rowvec elaps_t_data = arma::rowvec(elaps_times) * 1000; // convert to ms
  if (!elaps_t_data.empty())
  {
    double mean_elaps_t = arma::mean(elaps_t_data);
    double std_elaps_t = arma::stddev(elaps_t_data);
    double max_elaps_t = arma::max(elaps_t_data);
    double min_elaps_t = arma::min(elaps_t_data);

    std::cerr << "======= Elapsed time (ms) ======\n";
    std::cerr << "std_range: [" << std::max(0.,mean_elaps_t-std_elaps_t) << " -  " << mean_elaps_t + std_elaps_t <<"]\n";
    std::cerr << "mean     : " << mean_elaps_t << " +/- " << std_elaps_t <<"\n";
    std::cerr << "min      : " << min_elaps_t << "\n";
    std::cerr << "max      : " << max_elaps_t << "\n";
    std::cerr << "total_time: " << t << " sec\n";
    std::cerr << "==============================\n";
  }

  return Result(true);
}


// ===============================================================

void OnlineAdaptController::setCurrentPoseAsPickTarget(bool enable)
{
  const_pick_target.p = robot->getTaskPosition();
  const_pick_target.Q = robot->getTaskOrientation();

  if (enable)
  {
    getPickTargetPose = [this](arma::vec *p, arma::vec *Q)
    {
      *p = const_pick_target.p;
      *Q = const_pick_target.Q;
      return 1;
    };
  }
  else getPickTargetPose = [this](arma::vec *P, arma::vec *Q){ return pick_target_reader->getPose(P, Q); };
  
}

arma::vec OnlineAdaptController::getTaskWrench(const arma::vec Fext_prev) const
{
  arma::vec F_msr = robot->getTaskWrench();
  arma::vec Fext = ctrl_params.a_f*F_msr + (1 - ctrl_params.a_f)*Fext_prev;
  applyEnableDoFs(Fext);
  return Fext;
}

Result OnlineAdaptController::updateRobot()
{
  if (!robot->isOk())
    return Result(false, "The robot is not ok:\n" + robot->getErrMsg() + "\nAborting execution...");
  robot->update();
  return Result(true);
}

Result OnlineAdaptController::checkVelLimits(const arma::vec &V_cmd) const
{
  double v_norm = arma::norm(V_cmd.subvec(0, 2));
  if (v_norm > ctrl_params.vel_lim)
  {
    std::string msg = std::to_string(v_norm) + " > " + std::to_string(ctrl_params.vel_lim);
    return Result(false, "Velocity limit exceeded: " + msg);
  }
  return Result(true);
}

// ===============================================================

ExecResultMsg OnlineAdaptController::savePickModel(const std::string &filename)
{
  auto &model = pick_model;
  std::string err_msg;
  bool success = model.save(filename, &err_msg);
  if (success) return ExecResultMsg(ExecResultMsg::INFO, "Saved model successfully!");
  else return ExecResultMsg(ExecResultMsg::ERROR, err_msg);
}

ExecResultMsg OnlineAdaptController::loadPickModel(const std::string &filename)
{
  auto &model = pick_model;
  std::string err_msg;
  bool success = model.load(filename, &err_msg);
  if (success) return ExecResultMsg(ExecResultMsg::INFO, "Loaded model successfully!");
  else return ExecResultMsg(ExecResultMsg::ERROR, err_msg);
}

ExecResultMsg OnlineAdaptController::saveLoggedData()
{
  std::string count = std::to_string(save_counter) + "_";

  std::string f1 = getDefaultPath() + "/" + count + "exec_pick_frw.bin";
  std::string err_msg1 = "SUCCESS";
  bool success_1 = saveLoggedDataHelper(f1, frw_logger, &err_msg1);
  err_msg1 = "Save pick-forward: " + err_msg1;

  std::string f2 = getDefaultPath() + "/" + count + "exec_pick_rev.bin";
  std::string err_msg2 = "SUCCESS";
  bool success_2 = saveLoggedDataHelper(f2, rev_logger, &err_msg2);
  err_msg2 = "Save pick-reverse: " + err_msg2;

  std::string f3 = getDefaultPath() + "/" + count + "exec_place_frw.bin";
  std::string err_msg3 = "SUCCESS";
  bool success_3 = saveLoggedDataHelper(f3, frw_place_logger, &err_msg3);
  err_msg3 = "Save place-forward: " + err_msg3;

  std::string f4 = getDefaultPath() + "/" + count + "exec_place_rev.bin";
  std::string err_msg4 = "SUCCESS";
  bool success_4 = saveLoggedDataHelper(f4, rev_place_logger, &err_msg4);
  err_msg4 = "Save place-reverse: " + err_msg4;
  
  bool success = success_1 && success_2 && success_3 && success_4;
  std::string err_msg = err_msg1 + "\n" + err_msg2 + "\n" + err_msg3 + "\n" + err_msg4;

  save_counter++;

  if (success) return ExecResultMsg(ExecResultMsg::INFO, "Saved logged data successfully!");
  else return ExecResultMsg(ExecResultMsg::ERROR, "Error on saving logged data:\n" + err_msg);
}

bool OnlineAdaptController::saveLoggedDataHelper(const std::string &filename, const std::vector<LogState> &logger, std::string *err_msg)
{

  try
  {
    int n_data = logger.size();
    if (n_data == 0)
    {
      *err_msg = "The data are empty!";
      // throw std::runtime_error("The data are empty!");
      return true;
    }

    std::map<std::string, arma::mat> data;

    for (auto var_n : LogState::var_names) data[var_n] = arma::mat(log(var_n).size(), n_data);

    for (int j=0; j<n_data; j++)
    {
      auto &log_j = logger[j];
      for (auto var_n : LogState::var_names) data[var_n].col(j) = log_j(var_n);
    }

    io_::FileIO fid(filename, io_::FileIO::out | io_::FileIO::trunc);
    for (auto var_n : LogState::var_names) fid.write(var_n + "_data", data[var_n]);
    fid.close();
  }
  catch (std::exception &e)
  {
    *err_msg = e.what();
    // std::cerr << "\33[1;31m" << *err_msg << "\33[0m\n" << std::endl;
    return false;
  }

  return true;
}

ExecResultMsg OnlineAdaptController::loadParams()
{ 
  std::string filename = default_data_path + "ctrl_params.yaml";
  try{
    YAML::Node nh = YAML::LoadFile(filename);

    // ---------------- Model params ------------------
    pick_model.loadParams(filename, "model_params");
    place_model.loadParams(filename, "model_params");

    // ---------------- Controller params ------------------
    YAML::Node ctrl_node;
    if ( !YAML::getParam(nh, "ctrl_params", ctrl_node) ) throw std::runtime_error("Failed to load param \"ctrl_params\"...");
    if ( !YAML::getParam(ctrl_node, "Mp", ctrl_params.Mp) ) throw std::runtime_error("Failed to load param \"ctrl_params.Mp\"...");
    if ( !YAML::getParam(ctrl_node, "Dp", ctrl_params.Dp) ) throw std::runtime_error("Failed to load param \"ctrl_params.Dp\"...");
    if ( !YAML::getParam(ctrl_node, "Kp", ctrl_params.Kp) ) throw std::runtime_error("Failed to load param \"ctrl_params.Kp\"...");

    if ( !YAML::getParam(ctrl_node, "Mo", ctrl_params.Mo) ) throw std::runtime_error("Failed to load param \"ctrl_params.Mo\"...");
    if ( !YAML::getParam(ctrl_node, "Do", ctrl_params.Do) ) throw std::runtime_error("Failed to load param \"ctrl_params.Do\"...");
    if ( !YAML::getParam(ctrl_node, "Ko", ctrl_params.Ko) ) throw std::runtime_error("Failed to load param \"ctrl_params.Ko\"...");

    if ( !YAML::getParam(ctrl_node, "a_f", ctrl_params.a_f) ) throw std::runtime_error("Failed to load param \"ctrl_params.a_f\"...");

    if ( !YAML::getParam(ctrl_node, "vel_lim", ctrl_params.vel_lim) ) throw std::runtime_error("Failed to load param \"ctrl_params.vel_lim\"...");

    YAML::Node enabled_dofs_node;
    if ( !YAML::getParam(ctrl_node, "enabled_dofs", enabled_dofs_node) ) throw std::runtime_error("Failed to load param \"ctrl_params.enabled_dofs\"...");

    if ( !YAML::getParam(enabled_dofs_node, "dofs", ctrl_params.enabled_dofs.dofs) ) throw std::runtime_error("Failed to load param \"ctrl_params.enabled_dofs.dofs\"...");
    std::string frame;
    if ( !YAML::getParam(enabled_dofs_node, "frame", frame) ) throw std::runtime_error("Failed to load param \"ctrl_params.enabled_dofs.frame\"...");
    ctrl_params.enabled_dofs.frame = frame;
    if (frame.compare("base") == 0) applyEnableDoFs = [this](arma::vec &Fext){ Fext = Fext % ctrl_params.enabled_dofs.dofs; };
    else if (frame.compare("tool") == 0)
    {
      applyEnableDoFs = [this](arma::vec &Fext)
      {
        arma::mat R = robot->getTaskRotMat();
        Fext.subvec(0,2) = R.t() * Fext.subvec(0,2);
        Fext.subvec(3,5) = R.t() * Fext.subvec(3,5);
        Fext = Fext % ctrl_params.enabled_dofs.dofs;
        Fext.subvec(0,2) = R * Fext.subvec(0,2);
        Fext.subvec(3,5) = R * Fext.subvec(3,5);
      };
    }
    else throw std::runtime_error("Unrecognized frame '" + frame + "' for enabled dofs");

    // ---------------- Visualization params ------------------
    viz->loadParams(filename, "viz_params");

    // ----------------------------------------------------

    if ( !YAML::getParam(nh, "pick_model_filename", pick_model_filename) ) throw std::runtime_error("Failed to load param \"pick_model_filename\"...");
    pick_model_filename = getDefaultPath() + "/" + pick_model_filename;
    if ( !YAML::getParam(nh, "place_model_filename", place_model_filename) ) throw std::runtime_error("Failed to load param \"place_model_filename\"...");
    place_model_filename = getDefaultPath() + "/" + place_model_filename;

    if ( !YAML::getParam(nh, "pick_Tf", pick_Tf) ) throw std::runtime_error("Failed to load param \"pick_Tf\"...");
    if ( !YAML::getParam(nh, "place_Tf", place_Tf) ) throw std::runtime_error("Failed to load param \"place_Tf\"...");

    // if ( !YAML::getParam(nh, "place_vp_offsets", place_vp_offsets) ) throw std::runtime_error("Failed to load param \"place_vp_offsets\"...");

    if ( !YAML::getParam(nh, "a_Tf", a_Tf) ) throw std::runtime_error("Failed to load param \"a_Tf\"...");
    if ( !YAML::getParam(nh, "vp_dist_thres", vp_dist_thres) ) throw std::runtime_error("Failed to load param \"vp_dist_thres\"...");
    if ( !YAML::getParam(nh, "log_reserve", log_reserve) ) throw std::runtime_error("Failed to load param \"log_reserve\"...");
    if ( !YAML::getParam(nh, "log_every", log_every) ) throw std::runtime_error("Failed to load param \"log_every\"...");
    if ( !YAML::getParam(nh, "auto_save", auto_save) ) throw std::runtime_error("Failed to load param \"auto_save\"...");

    if ( !YAML::getParam(nh, "skip_pick", skip_pick) ) throw std::runtime_error("Failed to load param \"skip_pick\"...");

    if ( !YAML::getParam(nh, "phase_stop_tag_id", phase_stop_tag_id) ) throw std::runtime_error("Failed to load param \"phase_stop_tag_id\"...");
    if ( !YAML::getParam(nh, "adapt_to_robot_tag_id", adapt_to_robot_tag_id) ) throw std::runtime_error("Failed to load param \"adapt_to_robot_tag_id\"...");
    if ( !YAML::getParam(nh, "stop_exec_tag_id", stop_exec_tag_id) ) throw std::runtime_error("Failed to load param \"stop_exec_tag_id\"...");

    if ( !YAML::getParam(nh, "P_park", P_park) ) throw std::runtime_error("Failed to load param \"P_park\"...");
    if ( !YAML::getParam(nh, "Q_park", Q_park) ) throw std::runtime_error("Failed to load param \"Q_park\"...");

    if ( !YAML::getParam(nh, "obst_bounds", obst_bounds) ) throw std::runtime_error("Failed to load param \"obst_bounds\"...");
    arma::vec obst_extra_bounds;
    if ( !YAML::getParam(nh, "obst_extra_bounds", obst_extra_bounds) ) throw std::runtime_error("Failed to load param \"obst_extra_bounds\"...");
    obst_bounds += obst_extra_bounds;

    // --------------- Tag params --------------------

    pick_target_reader->loadParams(filename, "pick_target_params");
    vp_reader->loadParams(filename, "viapoint_params");
    place_target_reader->loadParams(filename, "place_target_params");
    obstacle_reader->loadParams(filename, "obstacle_params");

    // ----------------------------------------------------

    return ExecResultMsg(ExecResultMsg::INFO, "Loaded params successfully!");
  }
  catch (std::exception &e)
  { return ExecResultMsg(ExecResultMsg::ERROR, std::string("Failed to load params:") + e.what()); }
  
}

void OnlineAdaptController::setVizualizePickTarget(bool set, unsigned pub_rate)
{
  viz->viewTag("target", set, [this](arma::vec *p, arma::vec *Q){ return getPickTargetPose(p, Q); }, pub_rate);
}

void OnlineAdaptController::setVizualizeViapoint(bool set, unsigned pub_rate)
{
  viz->viewTag("via-point", set, [this](arma::vec *p, arma::vec *Q){ return vp_reader->getPoseHelper(p, Q); }, pub_rate);
}

void OnlineAdaptController::setVizualizePlaceTarget(bool set, unsigned pub_rate)
{
  viz->viewTag("place-target", set, [this](arma::vec *p, arma::vec *Q){ return place_target_reader->getPoseHelper(p, Q); }, pub_rate);
}

void OnlineAdaptController::setVizualizeObst(bool set, unsigned pub_rate)
{
  viz->viewTag("obst", set, [this](arma::vec *p, arma::vec *Q){ return obstacle_reader->getPoseHelper(p, Q); }, pub_rate);
}

void OnlineAdaptController::showObstacle(bool set)
{
  if (show_obst.read() == set) return;

  show_obst.set(set);
  if (!set) return;
  
  std::thread([this]()
  {
    arma::vec P_obst, Q_obst;
    while (show_obst())
    {
      if (obstacle_reader->getPose(&P_obst, &Q_obst)) viz->showObstacle(P_obst, Q_obst, 2*obst_bounds(0), 2*obst_bounds(1), 2*obst_bounds(2));
      std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
    viz->clearObstacle();
  }).detach();
}


bool OnlineAdaptController::reloadAprilTagListener(std::string *err_msg)
{
  std::string msg;
  if (!initAprilTagListenerFromFile(default_data_path + "/apriltag_listener.yaml", &msg))
  {
    if (err_msg) *err_msg = msg;
    if (gui) emit gui->showMsgSignal(ExecResultMsg(ExecResultMsg::ERROR, msg));
    return false;
  }
  return true;
}

bool OnlineAdaptController::initAprilTagListenerFromFile(const std::string &tag_listener_cfg_file, std::string *err_msg)
{
  arma::mat Tf_robot_cam;
  std::string Tf_parent_link;
  std::string Tf_child_link;
  std::string camera_base_link;
  bool publish_detections_to_tf;

  try
  {
    YAML::Node config = YAML::LoadFile(tag_listener_cfg_file);

    if ( !YAML::getParam(config, "publish_detections_to_tf", publish_detections_to_tf) )
      throw std::runtime_error(OnlineAdaptController_fun_+"Failed to load param 'publish_detections_to_tf'...");

    // ------------  parse 'camera_tf'  --------------
    YAML::Node camera_tf_node;
    if ( !YAML::getParam(config, "camera_tf", camera_tf_node) )
      throw std::runtime_error(OnlineAdaptController_fun_+"Failed to load param 'camera_tf'...");
    if ( !camera_tf_node.IsMap() )
        throw std::runtime_error(OnlineAdaptController_fun_+"'camera_tf' must be a struct...");

    YAML::Node tf_node;
    arma::vec tf_pos;
    arma::vec tf_quat;
    if ( !YAML::getParam(camera_tf_node, "transform", tf_node) )
      throw std::runtime_error(OnlineAdaptController_fun_+"Failed to load param 'camera_tf.transform'...");
    if ( !YAML::getParam(tf_node, "pos", tf_pos) )
      throw std::runtime_error(OnlineAdaptController_fun_+"Failed to load param 'camera_tf.transform.pos'...");
    if ( !YAML::getParam(tf_node, "orient", tf_quat) )
      throw std::runtime_error(OnlineAdaptController_fun_+"Failed to load param 'camera_tf.transform.orient'...");
    Tf_robot_cam = arma::mat().eye(4, 4);
    Tf_robot_cam.submat(0, 0, 2, 2) = math_::quat2rotm(tf_quat);
    Tf_robot_cam.submat(0, 3, 2, 3) = tf_pos;

    if ( !YAML::getParam(camera_tf_node, "parent_link", Tf_parent_link) )
      online_adapt_::PRINT_WARNING_MSG(OnlineAdaptController_fun_+"Failed to load param 'camera_tf.parent_link'...");

    if ( !YAML::getParam(camera_tf_node, "child_link", Tf_child_link) )
      online_adapt_::PRINT_WARNING_MSG(OnlineAdaptController_fun_+"Failed to load param 'camera_tf.child_link'...");
    if ( !YAML::getParam(camera_tf_node, "camera_base_link", camera_base_link) )
      online_adapt_::PRINT_WARNING_MSG(OnlineAdaptController_fun_+"Failed to load param 'camera_tf.camera_base_link'...");
  }
  catch (std::exception &e)
  {
    if (err_msg) *err_msg = OnlineAdaptController_fun_ + "[Apriltag listener will not be created]: " + e.what();
    return false;
  }

  // updateCameraPoseInRviz(Tf_robot_cam, Tf_parent_link, camera_base_link);

  tag_listener.reset( new apriltag_ros::AprilTagListener(tag_listener_cfg_file, "operations_tags_map") );
  tag_listener->setTagsTrasform(Tf_robot_cam, Tf_parent_link, Tf_child_link);
  tag_listener->publishDetectionsToTf(publish_detections_to_tf);

  return true;
}

void OnlineAdaptController::updateCameraPoseInRviz(const arma::mat &Tf_base_camOpt, const std::string &base, const std::string &cam_base)
{
  tf::Transform T_base_camOpt;
  T_base_camOpt.setOrigin( tf::Vector3( Tf_base_camOpt(0,3), Tf_base_camOpt(1,3), Tf_base_camOpt(2,3) ) );
  T_base_camOpt.setBasis( tf::Matrix3x3( Tf_base_camOpt(0,0), Tf_base_camOpt(0,1), Tf_base_camOpt(0,2), \
                                         Tf_base_camOpt(1,0), Tf_base_camOpt(1,1), Tf_base_camOpt(1,2), \
                                         Tf_base_camOpt(2,0), Tf_base_camOpt(2,1), Tf_base_camOpt(2,2) ) );
  
  // from urdf of realsense
  tf::StampedTransform T_camBase_camOpt;
  T_camBase_camOpt.setOrigin( tf::Vector3(0.0106, 0.0325, 0.0125) );
  T_camBase_camOpt.setRotation( tf::Quaternion(-0.5, 0.5, -0.5, 0.5) );

  tf::Transform T_base_camBase = T_base_camOpt * T_camBase_camOpt.inverse();
  tf::TransformBroadcaster().sendTransform(tf::StampedTransform(T_base_camBase, ros::Time::now(), base, cam_base));
}