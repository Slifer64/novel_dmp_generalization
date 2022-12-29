#ifndef AS64_gmp_lib_UTILS_H
#define AS64_gmp_lib_UTILS_H

//#define GMP_DEBUG_
//#define WSoG_DEBUG_
//#define GMP_DEBUG_
//#define GMPo_DEBUG_

namespace as64_
{

  namespace gmp_
  {

    enum UPDATE_TYPE
    {
      POS = 0,
      VEL = 1,
      ACCEL = 2
    };

    struct Phase
    {
      double x; ///< phase variable
      double x_dot; ///< phase variable 1st time derivative
      double x_ddot; ///< phase variable 2nd time derivative

      /** constructor.
       * @param[in] x: phase variable.
       * @param[in] x_dot: phase variable 1st time derivative.
       * @param[in] x_ddot: phase variable 2nd time derivative.
       */
      Phase(double x1, double x1_dot=0, double x1_ddot=0): x(x1), x_dot(x1_dot), x_ddot(x1_ddot) {}
    };

  } // namespace gmp_

} // namespace as64_

#endif // AS64_gmp_lib_UTILS_H
