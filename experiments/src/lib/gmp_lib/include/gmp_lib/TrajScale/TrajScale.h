#ifndef GMP_TRAJ_SCALE_H
#define GMP_TRAJ_SCALE_H

// Trajectory scaling class
// Calculates the scaling matrix for scaling spatially a trajectory from a
// nominal start/final position {Y0d, Ygd} to a new start/final position {Y0, Yg}.
// The scaling is performed according to the scaling method (see @TrajScale::ScaleMethod enum)
// The scaling is returned as matrix T_sc that premultiplies a position to
// scale it, i.e. Y_scaled = T_sc * Y
// The scaling is also invertible, i.e. Y = inv(T_sc) * Y_scaled

#include <memory>
#include <armadillo>

namespace as64_
{

namespace gmp_
{

class TrajScale
{
// ===================================
// =======  Public Functions  ========
// ===================================
public:

  typedef std::shared_ptr<TrajScale> Ptr;

  enum ScaleType
  {
    NONE = -1,
    PROP_SCALE = 0,
    ROT_MIN_SCALE = 1,
    ROT_WB_SCALE = 2
  };

  /** Constructor.
   * @param[in] n_dof: degrees of freedom.
   */
  TrajScale(unsigned n_dof);

  /** Returns the matrix that scales to the new init/target position, according to the set scaling method.
   * @return: The scaling matrix to the new init/target position.
   */
  arma::mat getScaling()
  {
    return this->T_sc;
  }

  /** Returns the inverse of the matrix that scales to the new init/target position, according to the set scaling method.
   * @return: The inverse scaling matrix to the new init/target position.
   */
  arma::mat getInvScaling()
  {
    return this->inv_T_sc;
  }

  /** Set the new start and final position.
   * @param[in] Y0: new start position.
   * @param[in] Yg: new final position.
   */
  void setNewStartFinalPos(const arma::vec &Y0, const arma::vec &Yg);

  /** Set the nominal start and final position.
   * @param[in] Y0d: nominal start position.
   * @param[in] Ygd: nominal final position.
   */
  void setNominalStartFinalPos(const arma::vec &Y0d, const arma::vec &Ygd);

  /** Returns the number of DoFs.
   * @return: number of DoFs.
   */
  unsigned numOfDoFs() const
  {
    return this->Y0d.size();
  }

  gmp_::TrajScale::Ptr deepCopy() const;

  // ==============================================
  // =======  Abstract (public) Functions  ========
  // ==============================================
  virtual enum ScaleType getScaleType() const = 0;

protected:

  // =================================================
  // =======  Abstract (protected) Functions  ========
  // =================================================

  virtual arma::mat calcScaling() const = 0;
  virtual arma::mat calcInvScaling() const = 0;

  // =======================================
  // =======  Protected Properties  ========
  // =======================================

  unsigned n_dof; ///< number of DoFs
  
  arma::vec Y0d; ///< Nominal start position.
  arma::vec Ygd; ///< Nominal final position.
  arma::vec Y0; ///< New start position.
  arma::vec Yg; ///< New final position.

  arma::mat T_sc; //scaling matrix, calculated based on the set scale method
  arma::mat inv_T_sc; //inverse scaling matrix, calculated based on the set scale method

}; // TrajScale

} // namespace gmp_

} // namespace as64_

#endif // GMP_TRAJ_SCALE_H
