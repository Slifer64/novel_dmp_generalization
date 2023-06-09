#ifndef $_PROJECT_384$_ONLINE_ADAPT_CONTROLLER_H
#define $_PROJECT_384$_ONLINE_ADAPT_CONTROLLER_H

#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>
#include <armadillo>

#include <online_adapt_controller/gui.h>


#include <robot_wrapper/robot.h>
#include <main_controller/main_controller.h>
#include <main_controller/controller.h>

#include <apriltag_ros/apriltag_listener.h>

#include <online_adapt_controller/utils/DMP_pp.h>
#include <online_adapt_controller/utils/visualizer.h>
#include <online_adapt_controller/utils/print.h>
#include <online_adapt_controller/utils/gripper.h>
#include <online_adapt_controller/utils/tag_reader.h>

using namespace as64_;

class LogState
{
public:

  LogState()
  {
    for (auto v_name : var_names) m[v_name] = arma::vec();
  }

  static std::vector<std::string> var_names;

  arma::vec &operator()(const std::string &var_name)
  {
    auto it = m.find(var_name);
    if (it == m.end()) throw std::runtime_error("[LogState::operator()]: Invalid member name \"" + var_name + "\"....\n");
    return it->second;
  }

  arma::vec operator()(const std::string &var_name) const
  {
    auto it = m.find(var_name);
    if (it == m.end()) throw std::runtime_error("[LogState::operator()]: Invalid member name \"" + var_name + "\"....\n");
    return it->second;
  }

private:

  std::map<std::string, arma::vec> m;
};

struct Result
{
  Result(bool success=true, const std::string &err_msg="")
  {
    this->success = success;
    this->err_msg = err_msg;
  }
  bool success;
  std::string err_msg;
};

class OnlineAdaptController : public Controller
{
public:

  OnlineAdaptController(MainController *main_ctrl, const std::string &ctrl_name);
  ~OnlineAdaptController();

  QPushButton *createGui(MainWindow *parent) override;

  std::string getDefaultPath() const { return default_data_path; }

  double getProgress() const { return x_progress; }

  void setVizualizePickTarget(bool set, unsigned pub_rate);
  void setVizualizeViapoint(bool set, unsigned pub_rate);
  void setVizualizePlaceTarget(bool set, unsigned pub_rate);
  void setVizualizeObst(bool set, unsigned pub_rate);

  ExecResultMsg saveLoggedData();

  bool saveLoggedDataHelper(const std::string &filename, const std::vector<LogState> &logger, std::string *err_msg);

  ExecResultMsg savePickModel(const std::string &filename);

  ExecResultMsg loadPickModel(const std::string &filename);

public:

  bool keep_alive=true;
  void launchExternalGripperCtrlThread();

  ros::NodeHandle node;

  OnlineAdaptWin *gui = 0;

  std::unique_ptr<apriltag_ros::AprilTagListener> tag_listener;
  bool reloadAprilTagListener(std::string *err_msg=0);
  bool initAprilTagListenerFromFile(const std::string &tag_listener_cfg_file, std::string *err_msg=0);
  void updateCameraPoseInRviz(const arma::mat &Tf_base_camOpt, const std::string &base, const std::string &cam_base);

  bool getTagPose(int tag_id, arma::vec *P, arma::vec *Q) const;

  std::unique_ptr<online_adapt_::TagReader> pick_target_reader;
  std::unique_ptr<online_adapt_::TagReader> vp_reader;
  std::unique_ptr<online_adapt_::TagReader> place_target_reader;
  std::unique_ptr<online_adapt_::TagReader> obstacle_reader;

  arma::vec obst_bounds;
  void showObstacle(bool set);
  thr_::MtxVar<bool> show_obst;

  void start();
  ExecResultMsg stop();

  Result init();

  Result executePick(bool reverse);

  Result executePlace(bool reverse);

  ExecResultMsg loadParams();

  void setCurrentPoseAsPickTarget(bool enable);
  struct
  {
    arma::vec p;
    arma::vec Q;
  } const_pick_target;
  

  std::function<int(arma::vec *, arma::vec *)> getPickTargetPose;

  bool run_;

  arma::vec getTaskWrench(const arma::vec Fext_prev) const;
  Result updateRobot();
  Result checkVelLimits(const arma::vec &V) const;

  arma::vec P0, Q0;

  arma::vec P_park;
  arma::vec Q_park;

  struct
  {
    arma::vec P;
    arma::vec Q;
  } current_pose;

  // ------- Model ---------

  DMP_pp pick_model;
  DMP_pp place_model;
  std::string pick_model_filename;
  std::string place_model_filename;

  std::unique_ptr<online_adapt_::Gripper> gripper;

  // -----------------------

  std::string default_data_path;

  thr_::Semaphore exec_stop_sem;

  std::function<void(arma::vec &)> applyEnableDoFs;

  struct{
    double Mp;
    double Dp;
    double Kp;

    double Mo;
    double Do;
    double Ko;

    double a_f; // filtering coeff for wrench (1: no filtering)

    struct{
      arma::uvec dofs;
      std::string frame;
    } enabled_dofs;

    double vel_lim;
  } ctrl_params;

  double pick_Tf;
  double place_Tf;
  double a_Tf;

  // std::vector<double> place_vp_offsets;

  thr_::MtxVar<bool> enable_model_adapt;
  thr_::MtxVar<bool> phase_stop;
  thr_::MtxVar<bool> reset_model;

  unsigned long log_reserve;
  std::vector<LogState> frw_logger;
  std::vector<LogState> rev_logger;
  std::vector<LogState> frw_place_logger;
  std::vector<LogState> rev_place_logger;
  int log_every;
  LogState log;
  bool auto_save;
  int save_counter=0;

  double vp_dist_thres;

  bool sim;

  std::string base_link;
  std::string marker_array_topic;

  
  // ------- Visualization ---------
  std::unique_ptr<online_adapt_::Visualizer> viz;

  // -------------------------------

  double x_progress;

  bool execute_reverse = false;
  bool execute_place = false;
  bool retract_from_place = false;

  bool skip_pick=false;

  int phase_stop_tag_id;
  int adapt_to_robot_tag_id;
  int stop_exec_tag_id;
  void launchTagCtrlSignalsThread();
  std::thread tag_ctrl_thead;
};

#endif // $_PROJECT_384$_ONLINE_ADAPT_CONTROLLER_H
