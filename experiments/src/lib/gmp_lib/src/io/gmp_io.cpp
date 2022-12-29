#include <gmp_lib/io/gmp_io.h>
#include <gmp_lib/GMP/GMP_Update.h>

namespace as64_
{

namespace gmp_
{

  void write(const gmp_::GMP *gmp, const std::string &filename, const std::string &prefix)
  {
    gmp_::FileIO fid(filename, gmp_::FileIO::out|gmp_::FileIO::trunc);
    write(gmp, fid, prefix);
  }

  void write(const gmp_::GMP *gmp, gmp_::FileIO &fid, const std::string &prefix)
  {
    unsigned N_kernels = gmp->numOfKernels();
    unsigned n_dofs = gmp->numOfDoFs();

    fid.write(prefix + "weights", gmp->W);
    fid.write(prefix + "damping", gmp->D);
    fid.write(prefix + "stiffness", gmp->K);
    fid.write(prefix + "N_kernels", N_kernels);
    fid.write(prefix + "N_DoFs", n_dofs);
    fid.write(prefix + "scale_type", (int)(gmp->traj_sc->getScaleType()));
    fid.write(prefix + "c", gmp->getCenters());
    fid.write(prefix + "h", gmp->getInvWidths());
  }

  void read(gmp_::GMP *gmp, const std::string &filename, const std::string &prefix)
  {
    gmp_::FileIO fid(filename, gmp_::FileIO::in);
    read(gmp, fid, prefix);
  }

  void read(gmp_::GMP *gmp, gmp_::FileIO &fid, const std::string &prefix)
  {
    fid.read(prefix + "weights", gmp->W);
    fid.read(prefix + "damping", gmp->D);
    fid.read(prefix + "stiffness", gmp->K);
    int scale_type;
    fid.read(prefix + "scale_type", scale_type);

    arma::vec c, h;
    fid.read(prefix + "c", c);
    fid.read(prefix + "h", h);
    gmp->setKernels(c, h);

    unsigned n_dofs = gmp->numOfDoFs();

    gmp->Y0d = gmp->W*gmp->regressVec(0);
    gmp->Ygd = gmp->W*gmp->regressVec(1);
    gmp->Y0 = gmp->Y0d;
    gmp->Yg = gmp->Ygd;
    gmp->y_dot = arma::vec().zeros(n_dofs);
    gmp->z_dot = arma::vec().zeros(n_dofs);

    if (scale_type == TrajScale::PROP_SCALE) gmp->setScaleMethod( gmp_::TrajScale::Ptr(new gmp_::TrajScale_Prop(n_dofs)) );
    else if (scale_type == TrajScale::ROT_MIN_SCALE) gmp->setScaleMethod( gmp_::TrajScale::Ptr(new gmp_::TrajScale_Rot_min()) );
    else if (scale_type == TrajScale::ROT_WB_SCALE) gmp->setScaleMethod( gmp_::TrajScale::Ptr(new gmp_::TrajScale_Rot_wb()) );
    else if (scale_type == TrajScale::NONE) gmp->setScaleMethod( gmp_::TrajScale::Ptr(new gmp_::TrajScale_None(n_dofs)) );
    else throw std::runtime_error("[read]: Unsupported scale type \"" + std::to_string(scale_type) + "\"...\n");
  }


  void write(const gmp_::GMPo *gmp, const std::string &filename, const std::string &prefix)
  {
    gmp_::FileIO fid(filename, gmp_::FileIO::out|gmp_::FileIO::trunc);
    write(gmp, fid, prefix);
  }

  void write(const gmp_::GMPo *gmp, gmp_::FileIO &fid, const std::string &prefix)
  {
    write(dynamic_cast<const GMP *>(gmp), fid, prefix);
    fid.write(prefix + "Qd0", gmp->Qd0);
  }

  void read(gmp_::GMPo *gmp, const std::string &filename, const std::string &prefix)
  {
    gmp_::FileIO fid(filename, gmp_::FileIO::in);
    read(gmp, fid, prefix);
  }

  void read(gmp_::GMPo *gmp, gmp_::FileIO &fid, const std::string &prefix)
  {
    read(dynamic_cast<GMP *>(gmp), fid, prefix);
    fid.read(prefix + "Qd0", gmp->Qd0);
    gmp->setQ0(gmp->Qd0);
  }

} // namespace gmp_

} // namespace as64_
