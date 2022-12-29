#include <osqp_lib/quadprog.h>

#include <thread>
#include <string>
#include <iostream>

namespace osqp_
{

template<typename T>
T *vector_deepCopy(const T *src, c_int n)
{
  // return const_cast<T *>(src);

  T *dest = (T *)c_malloc(n*sizeof(T));
  memcpy( (void *)dest, (void *)src, sizeof(T)*n);
  return dest;
}

csc *csc_deepCopy(const csc *obj)
{
  return copy_csc_mat(obj);
}

OSQPPolish *OSQPPolish_deepCopy(const OSQPPolish *obj, c_int n, c_int m)
{
  // return const_cast<OSQPPolish *>(obj);
  
  OSQPPolish *cp = (OSQPPolish *)c_malloc(sizeof(OSQPPolish));

  // make shallow copy first
  // *cp = *obj;

  cp->n_low = obj->n_low;
  cp->n_upp = obj->n_upp;

  cp->obj_val = obj->obj_val;
  cp->pri_res = obj->pri_res;
  cp->dua_res = obj->dua_res;

  // then make deepCopy for the required attributes

  // cp->Ared = csc_deepCopy(obj->Ared);

  cp->A_to_Alow = vector_deepCopy(obj->A_to_Alow, m);
  cp->A_to_Aupp = vector_deepCopy(obj->A_to_Aupp, m);
  cp->Alow_to_A = vector_deepCopy(obj->Alow_to_A, m);
  cp->Aupp_to_A = vector_deepCopy(obj->Aupp_to_A, m);

  cp->x = vector_deepCopy(obj->x, n);
  cp->y = vector_deepCopy(obj->y, m);
  cp->z = vector_deepCopy(obj->z, m);

  return cp;
}


LinSysSolver *LinSysSolver_deepCopy(const LinSysSolver *obj)
{
  return const_cast<LinSysSolver *>(obj);

  LinSysSolver *cp = (LinSysSolver *)c_malloc(sizeof(LinSysSolver));
  *cp = *obj;

  return cp;
}

OSQPTimer *OSQPTimer_deepCopy(const OSQPTimer *obj)
{
  // return const_cast<OSQPTimer *>(obj);

  OSQPTimer *cp = (OSQPTimer *)c_malloc(sizeof(OSQPTimer));

  *cp = *obj;

  return cp;
}

OSQPData *OSQPData_deepCopy(const OSQPData *obj)
{
  // return const_cast<OSQPData *>(obj);

  OSQPData *cp = (OSQPData *)c_malloc(sizeof(OSQPData));

  cp->n = obj->n; // number of variables n
  cp->m = obj->m; // number of constraints m
  cp->P = csc_deepCopy(obj->P); // the upper triangular part of the quadratic cost matrix P in csc format (size n x n).
  cp->A = csc_deepCopy(obj->A); // linear constraints matrix A in csc format (size m x n)

  cp->q = vector_deepCopy(obj->q, obj->n); // dense array for linear part of cost function (size n)

  cp->l = vector_deepCopy(obj->l, obj->m); // dense array for lower bound (size m)
  cp->u = vector_deepCopy(obj->u, obj->m); // dense array for upper bound (size m)

  return cp;
}

OSQPSettings *OSQPSettings_deepCopy(const OSQPSettings *obj)
{
  //return const_cast<OSQPSettings *>(obj);

  // return copy_settings(obj);

  OSQPSettings *cp = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
  *cp = *obj;

  return cp;
}

OSQPScaling *OSQPScaling_deepCopy(const OSQPScaling *obj, c_int n, c_int m)
{
  // return const_cast<OSQPScaling *>(obj);

  OSQPScaling *cp = (OSQPScaling *)c_malloc(sizeof(OSQPScaling));

  cp->c = obj->c;
  cp->D = vector_deepCopy(obj->D, n);
  cp->E = vector_deepCopy(obj->E, m);

  cp->cinv = obj->cinv;
  cp->Dinv = vector_deepCopy(obj->Dinv, n);
  cp->Einv = vector_deepCopy(obj->Einv, m);

  return cp;
}


OSQPSolution *OSQPSolution_deepCopy(const OSQPSolution *obj, c_int n, c_int m)
{
  // return const_cast<OSQPSolution *>(obj);

  OSQPSolution *cp = (OSQPSolution *)c_malloc(sizeof(OSQPSolution));

  cp->x = vector_deepCopy(obj->x, n); // primal solution
  cp->y = vector_deepCopy(obj->y, m); // Lagrange multiplier associated to \f$l <= Ax <= u\f$

  return cp;
}


OSQPInfo *OSQPInfo_deepCopy(const OSQPInfo *obj)
{
  // return const_cast<OSQPInfo *>(obj);

  OSQPInfo *cp = (OSQPInfo *)c_malloc(sizeof(OSQPInfo));
  *cp = *obj;

  return cp;
}


OSQPWorkspace *OSQPWorkspace_deepCopy(const OSQPWorkspace *work)
{
  // return const_cast<OSQPWorkspace *>(work);

  OSQPWorkspace *cp = (OSQPWorkspace *)c_malloc(sizeof(OSQPWorkspace));

  // *cp = *work;

  cp->data = OSQPData_deepCopy(work->data);

  c_int n = cp->data->n;
  c_int m = cp->data->m;

  // cp->linsys_solver = LinSysSolver_deepCopy(work->linsys_solver);
  cp->linsys_solver = work->linsys_solver;

  # ifndef EMBEDDED
    if (work->settings->polish) cp->pol = OSQPPolish_deepCopy(work->pol, n, m);
  # endif // ifndef EMBEDDED

  // Vector used to store a vectorized rho parameter
  cp->rho_vec = vector_deepCopy(work->rho_vec, m); // vector of rho values
  cp->rho_inv_vec = vector_deepCopy(work->rho_inv_vec, m); // vector of inv rho values

  # if EMBEDDED != 1
    cp->constr_type = vector_deepCopy(work->constr_type, m); // Type of constraints: loose (-1), equality (1), inequality (0)
  # endif // if EMBEDDED != 1

  // Iterates
  cp->x = vector_deepCopy(work->x, n); // Iterate x
  cp->y = vector_deepCopy(work->x, m); // Iterate y
  cp->z = vector_deepCopy(work->x, m); // Iterate z

  cp->xz_tilde = vector_deepCopy(work->xz_tilde, n+m); // Iterate xz_tilde
  
  cp->x_prev = vector_deepCopy(work->x_prev, n); // Previous x

  cp->z_prev = vector_deepCopy(work->z_prev, m); // Previous z

  // Primal and dual residuals workspace variables
  cp->Ax = vector_deepCopy(work->Ax, m); // scaled A * x
  cp->Px = vector_deepCopy(work->Px, n); // scaled P * x
  cp->Aty = vector_deepCopy(work->Aty, n); // scaled A' * y

  //Primal infeasibility variables
  cp->delta_y = vector_deepCopy(work->delta_y, m); // difference between consecutive dual iterates
  cp->Atdelta_y = vector_deepCopy(work->Atdelta_y, n); // A' * delta_y

  // Dual infeasibility variables
  cp->delta_x = vector_deepCopy(work->delta_x, n); // difference between consecutive primal iterates
  cp->Pdelta_x = vector_deepCopy(work->Pdelta_x, n); // P * delta_x
  cp->Adelta_x = vector_deepCopy(work->Adelta_x, m); // A * delta_x

  //Temporary vectors used in scaling
  cp->D_temp = vector_deepCopy(work->D_temp, n); // temporary primal variable scaling vectors
  cp->D_temp_A = vector_deepCopy(work->D_temp_A, n); // temporary primal variable scaling vectors storing norms of A columns
  cp->E_temp = vector_deepCopy(work->E_temp, m); // temporary constraints scaling vectors storing norms of A' columns

  cp->settings = OSQPSettings_deepCopy(work->settings); // problem settings
  cp->scaling = OSQPScaling_deepCopy(work->scaling, n, m); // scaling vectors
  cp->solution = OSQPSolution_deepCopy(work->solution, n, m); // problem solution
  cp->info = OSQPInfo_deepCopy(work->info); // solver information

  # ifdef PROFILING
    cp->timer = OSQPTimer_deepCopy(work->timer);
    cp->first_run = work->first_run;
    cp->clear_update_time = work->clear_update_time;
    cp->first_run = work->first_run;
    cp->rho_update_from_solve = work->rho_update_from_solve;
  # endif // ifdef PROFILING

  # ifdef PRINTING
    cp->summary_printed = work->summary_printed;
  # endif // ifdef PRINTING

  return cp;
}


void printQSQPSettings(const OSQPSettings *settings)
{
  std::cout << "###################################\n"
            << "########  QSQP - Settings  ########\n"
            << "###################################\n";

  std::cout << "rho: " << settings->rho << "\n"
            << "sigma: " << settings->sigma << "\n"
            << "scaling: " << settings->scaling << "\n";

  # if EMBEDDED != 1
    std::cout << "adaptive_rho:" << settings->adaptive_rho << "\n"
              << "adaptive_rho_interval:" << settings->adaptive_rho_interval << "\n"
              << "adaptive_rho_tolerance:" << settings->adaptive_rho_tolerance << "\n";
    #  ifdef PROFILING
      std::cout << "adaptive_rho_fraction:" << settings->adaptive_rho_fraction << "\n";
    #  endif // Profiling
  # endif // EMBEDDED != 1

  std::cout << "max_iter:" << settings->max_iter << "\n"
            << "eps_abs:" << settings->eps_abs << "\n"
            << "eps_rel:" << settings->eps_rel << "\n"
            << "eps_prim_inf:" << settings->eps_prim_inf << "\n"
            << "eps_dual_inf:" << settings->eps_dual_inf << "\n"
            << "alpha:" << settings->alpha << "\n";

  std::string lin_sys_solver_name;
  if (settings->linsys_solver == QDLDL_SOLVER) lin_sys_solver_name = "QDLDL_SOLVER";
  else if (settings->linsys_solver == MKL_PARDISO_SOLVER) lin_sys_solver_name = "MKL_PARDISO_SOLVER";
  std::cout << "linsys_solver:" << lin_sys_solver_name << "\n";

  # ifndef EMBEDDED
    std::cout << "delta:" << settings->delta << "\n"
              << "polish:" << settings->polish << "\n"
              << "polish_refine_iter:" << settings->polish_refine_iter << "\n"
              << "verbose:" << settings->verbose << "\n";
  # endif // ifndef EMBEDDED

  std::cout << "scaled_termination:" << settings->scaled_termination << "\n"
            << "check_termination:" << settings->check_termination << "\n"
            << "warm_start:" << settings->warm_start << "\n";

  # ifdef PROFILING
    std::cout << "time_limit:" << settings->time_limit << "\n";
  # endif // ifdef PROFILING
}

void printQSQPInfo(const OSQPWorkspace *work)
{
  // OSQPPolish *pol = work->pol;
  // std::cout << "#########  OSQPPolish  #########\n";
  // std::cout << "Ared: " << pol->Ared << "\n"
  //           << "n_low: " << pol->n_low << "\n"
  //           << "n_upp: " << pol->n_upp << "\n"
  //           << "A_to_Alow: " << pol->A_to_Alow << "\n"
  //           << "A_to_Aupp: " << pol->A_to_Aupp << "\n"
  //           << "Alow_to_A: " << pol->Alow_to_A << "\n"
  //           << "Aupp_to_A: " << pol->Aupp_to_A << "\n"
  //           << "x: " << pol->x << "\n"
  //           << "y: " << pol->y << "\n"
  //           << "z: " << pol->z << "\n"
  //           << "obj_val: " << pol->obj_val << "\n"
  //           << "pri_res: " << pol->pri_res << "\n"
  //           << "dua_res: " << pol->dua_res << "\n";
  // std::cout << "################################\n";

  // LinSysSolver *solver = work->linsys_solver;
  // std::cout << "#########  LinSysSolver  #########\n";
  // std::cout << "type: " << solver->type << "\n"
  //           << "solve: " << solver->solve << "\n"
  //           << "free: " << solver->free << "\n"
  //           << "update_matrices: " << solver->update_matrices << "\n"
  //           << "update_rho_vec: " << solver->update_rho_vec << "\n"
  //           << "nthreads: " << solver->nthreads << "\n";
  // std::cout << "################################\n";
}


QuadProgSolution quadprog(const arma::mat &H, const arma::mat &f, const arma::mat &A, const arma::mat &lb, const arma::mat &ub,
              const arma::mat &Aeq, const arma::mat &beq, const QuadProgOptions &options)
{
  arma::mat f_ = f;
  arma::mat lb_ = arma::join_vert(lb, beq);
  arma::mat ub_ = arma::join_vert(ub, beq);

  int n_var = A.n_cols;
  int n_constr = A.n_rows + Aeq.n_rows;

  CSC_mat P_cs(H, true);
  CSC_mat A_cs(arma::join_vert(A, Aeq));

  // Exitflag
  c_int exitflag = 0;

  // Workspace structures
  OSQPWorkspace *work;
  OSQPSettings  *settings = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
  OSQPData      *data     = (OSQPData *)c_malloc(sizeof(OSQPData));

  // Populate data
  if (data)
  {
    data->n = n_var;
    data->m = n_constr;
    data->P = csc_matrix(data->n, data->n, P_cs.nnz, P_cs.dataPtr(), P_cs.rowIndPtr(), P_cs.csPtr());
    data->A = csc_matrix(data->m, data->n, A_cs.nnz, A_cs.dataPtr(), A_cs.rowIndPtr(), A_cs.csPtr());
  }

  // Define solver settings as default
  if (settings)
  {
    osqp_set_default_settings(settings);
    // printQSQPSettings(settings);
    settings->alpha = 1.0; // Change alpha parameter
    settings->warm_start = options.warm_start;
    settings->polish = options.polish;
    settings->time_limit = options.time_limit;
    settings->max_iter = options.max_iters;
    settings->verbose = options.verbose;
  }

  int n_dim = f.n_cols;
  arma::mat W(n_var, n_dim);

  data->q = &(f_.col(0)(0));
  data->l = &(lb_.col(0)(0));
  data->u = &(ub_.col(0)(0));

  // ========  Setup workspace  =========
  OSQPWorkspace *work_temp;
  // exitflag = osqp_setup(&work_temp, data, settings);
  // if (exitflag != 0) throw std::runtime_error(std::string("[osqp_setup]: ") + work_temp->info->status + "\n");

  QuadProgSolution solution;
  solution.x.resize(n_var, n_dim);
  solution.status = QuadProgSolutionStatus::OPTIMAL;
  solution.exit_msg = "";

  // std::vector<OSQPWorkspace *> work_(n_dim);
  // work_[0] = work_temp;
  // for (int i=1; i<n_dim; i++)
  // {
  //   work_[i] = OSQPWorkspace_deepCopy(work_temp);

  //   exitflag = osqp_update_lin_cost(work_[i], &(f_.col(i)(0)));
  //   if (exitflag != 0) throw std::runtime_error(std::string("[osqp_setup]: ") + work_[i]->info->status + "\n");

  //   exitflag = osqp_update_bounds(work_[i], &(lb_.col(i)(0)), &(ub_.col(i)(0)));
  //   if (exitflag != 0) throw std::runtime_error(std::string("[osqp_setup]: ") + work_[i]->info->status + "\n");
  // }

  if (options.parallel)
  {

    std::vector<OSQPWorkspace *> work_(n_dim);
    for (int i=0; i<n_dim; i++)
    {
      data->q = &(f_.col(i)(0));
      data->l = &(lb_.col(i)(0));
      data->u = &(ub_.col(i)(0));

      exitflag = osqp_setup(&work_temp, data, settings);
      if (exitflag != 0) throw std::runtime_error(std::string("[osqp_setup]: ") + work_temp->info->status + "\n");

      work_[i] = work_temp;
    }


    std::vector<std::string> exit_msg(n_dim);
    std::vector<arma::vec> x_solution(n_dim);
    std::vector<QuadProgSolutionStatus> solution_status(n_dim);
    std::vector<std::thread> thr(n_dim);
    for (int k=0; k<n_dim; k++)
    {
      OSQPWorkspace *work = work_[k];
      thr[k] = std::thread( [work, k, &x_solution, &solution_status, &exit_msg]()
      {
        // Solve Problem
        osqp_solve(work);

        c_int n_var = work->data->n;

        c_int sol_status = work->info->status_val;
        if (sol_status == OSQP_SOLVED)
        {
          exit_msg[k] = std::string("Dim ") + std::to_string(k+1) + ": Function converged to the solution.";
          solution_status[k] = QuadProgSolutionStatus::OPTIMAL;
        }
        else if (sol_status == OSQP_SOLVED_INACCURATE || sol_status == OSQP_MAX_ITER_REACHED || sol_status == OSQP_TIME_LIMIT_REACHED)
        {
          exit_msg[k] = std::string("Dim ") + std::to_string(k+1) + ": " + work->info->status;
          solution_status[k] = QuadProgSolutionStatus::SUBOPTIMAL;
        }
        else
        {
          exit_msg[k] = std::string("Dim ") + std::to_string(k+1) + ": " + work->info->status;
          solution_status[k] = QuadProgSolutionStatus::FAILED;
        }
        x_solution[k] = arma::mat(work->solution->x, n_var, 1, true);

        osqp_cleanup(work);
      });

    }
    for (int i=0; i<n_dim; i++) thr[i].join();

    for (int i=0; i<n_dim; i++)
    {
      solution.x.col(i) = x_solution[i];
      solution.exit_msg += exit_msg[i] + "\n";
      if (solution_status[i] == QuadProgSolutionStatus::FAILED) 
        solution.status = QuadProgSolutionStatus::FAILED;
      if (solution.status != QuadProgSolutionStatus::FAILED && solution_status[i] == QuadProgSolutionStatus::SUBOPTIMAL)
        solution.status = QuadProgSolutionStatus::SUBOPTIMAL;
    }

  }
  else
  {
    
  }
 

/*
  for (int k=0; k<n_dim; k++)
  {
    OSQPWorkspace *work = work_[k];

    // Solve Problem
    osqp_solve(work);

    c_int sol_status = work->info->status_val;
    if (!solution.exit_msg.empty()) solution.exit_msg += "\n";
    if (sol_status == OSQP_SOLVED)
    {
      solution.exit_msg += std::string("Dim ") + std::to_string(k+1) + ": Function converged to the solution.";
      solution.x.col(k) = arma::mat(work->solution->x, n_var, 1, true);
    }
    else if (sol_status == OSQP_SOLVED_INACCURATE || sol_status == OSQP_MAX_ITER_REACHED || sol_status == OSQP_TIME_LIMIT_REACHED)
    {
      solution.exit_msg += std::string("Dim ") + std::to_string(k+1) + ": " + work->info->status;
      solution.status = QuadProgSolutionStatus::SUBOPTIMAL;
      solution.x.col(k) = arma::mat(work->solution->x, n_var, 1, true);
    }
    else
    {
      solution.exit_msg += std::string("Dim ") + std::to_string(k+1) + ": " + work->info->status;
      solution.status = QuadProgSolutionStatus::FAILED;
      // do not assign failed solution
    }
  }

  // for (int i=1; i<n_dim; i++) osqp_cleanup(work_[i]);
*/

/*
  work = work_temp;
  for (int k=0; k<n_dim; k++)
  {
    // Solve Problem
    osqp_solve(work);

    c_int sol_status = work->info->status_val;
    if (!solution.exit_msg.empty()) solution.exit_msg += "\n";
    if (sol_status == OSQP_SOLVED)
    {
      solution.exit_msg += std::string("Dim ") + std::to_string(k+1) + ": Function converged to the solution.";
      solution.x.col(k) = arma::mat(work->solution->x, n_var, 1, true);
    }
    else if (sol_status == OSQP_SOLVED_INACCURATE || sol_status == OSQP_MAX_ITER_REACHED || sol_status == OSQP_TIME_LIMIT_REACHED)
    {
      solution.exit_msg += std::string("Dim ") + std::to_string(k+1) + ": " + work->info->status;
      solution.status = QuadProgSolutionStatus::SUBOPTIMAL;
      solution.x.col(k) = arma::mat(work->solution->x, n_var, 1, true);
    }
    else
    {
      solution.exit_msg += std::string("Dim ") + std::to_string(k+1) + ": " + work->info->status;
      solution.status = QuadProgSolutionStatus::FAILED;
      // do not assign failed solution
    }

    if (k != n_dim-1)
    {
      osqp_update_lin_cost(work, &(f_.col(k + 1)(0)));
      osqp_update_bounds(work, &(lb_.col(k + 1)(0)), &(ub_.col(k + 1)(0)));
    }
  }
  osqp_cleanup(work);
*/


  // Cleanup
  if (data)
  {
    if (data->A) c_free(data->A);
    if (data->P) c_free(data->P);
    c_free(data);
  }
  if (settings) c_free(settings);

  return solution;
}

} // namespace osqp_
