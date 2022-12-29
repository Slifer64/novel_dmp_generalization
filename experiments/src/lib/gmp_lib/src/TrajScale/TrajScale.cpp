#include <gmp_lib/TrajScale/TrajScale.h>

#include <gmp_lib/TrajScale/TrajScale_Prop.h>
#include <gmp_lib/TrajScale/TrajScale_Rot_min.h>
#include <gmp_lib/TrajScale/TrajScale_Rot_wb.h>
#include <gmp_lib/TrajScale/TrajScale_None.h>

namespace as64_
{

namespace gmp_
{


TrajScale::TrajScale(unsigned n_dof)
{
  this->n_dof = n_dof;

  this->Y0d = arma::vec().zeros(n_dof);
  this->Ygd = arma::vec().ones(n_dof);
  //Initialization required for 'setScaleMethod'
  //will be overwirtten then on the 'setInitFinalPos' call
  this->Y0 = arma::vec().zeros(n_dof);
  this->Yg = arma::vec().ones(n_dof);
}


void TrajScale::setNewStartFinalPos(const arma::vec &Y0, const arma::vec &Yg)
{
  this->Y0 = Y0;
  this->Yg = Yg;

  this->T_sc = this->calcScaling();
  this->inv_T_sc = this->calcInvScaling();
}

void TrajScale::setNominalStartFinalPos(const arma::vec &Y0d, const arma::vec &Ygd)
{
  this->Y0d = Y0d;
  this->Ygd = Ygd;

  this->T_sc = this->calcScaling();
  this->inv_T_sc = this->calcInvScaling();
}

gmp_::TrajScale::Ptr TrajScale::deepCopy() const
{
  gmp_::TrajScale::Ptr cp_obj;

  arma::vec n_wb;

  gmp_::TrajScale::ScaleType sc_t = this->getScaleType();
  switch (getScaleType())
  {
    case gmp_::TrajScale::PROP_SCALE:
      cp_obj.reset( new gmp_::TrajScale_Prop(numOfDoFs()) );
      break;
    case gmp_::TrajScale::ROT_MIN_SCALE:
      cp_obj.reset( new gmp_::TrajScale_Rot_min() );
      break;
    case gmp_::TrajScale::ROT_WB_SCALE:
      cp_obj.reset( new gmp_::TrajScale_Rot_wb() );
      n_wb = dynamic_cast<const gmp_::TrajScale_Rot_wb *>(this)->getWorkBenchNormal();
      dynamic_cast<gmp_::TrajScale_Rot_wb *>(cp_obj.get())->setWorkBenchNormal(n_wb);
      break;
    case gmp_::TrajScale::NONE:
      cp_obj.reset( new gmp_::TrajScale_None(numOfDoFs()) );
      break;
  }

  cp_obj->setNominalStartFinalPos(Y0d, Ygd);
  cp_obj->setNewStartFinalPos(Y0, Yg);

  return cp_obj;
}

} // namespace gmp_

} // namespace as64_
