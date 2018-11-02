#include "cSMHdCKM_cSMHdCKMRHN_two_scale_matching.cpp"

namespace flexiblesusy {

void cSMHdCKM_cSMHdCKMRHN_matching::match_low_to_high_scale_model()
{
   upper->set_Yd(lower->get_Yd());
   upper->set_Ye(lower->get_Ye());
   upper->set_Yu(lower->get_Yu());
   upper->set_g1(lower->get_g1());
   upper->set_g2(lower->get_g2());
   upper->set_g3(lower->get_g3());
   upper->set_vd(lower->get_vd());
   upper->set_v(lower->get_v());
   upper->set_Lambdax(lower->get_Lambdax());
   upper->set_mu2(lower->get_mu2());

   upper->set_scale(lower->get_scale());
}

void cSMHdCKM_cSMHdCKMRHN_matching::match_high_to_low_scale_model()
{
   update_scale();

   lower->set_Yd(upper->get_Yd());
   lower->set_Ye(upper->get_Ye());
   lower->set_Yu(upper->get_Yu());
   lower->set_g1(upper->get_g1());
   lower->set_g2(upper->get_g2());
   lower->set_g3(upper->get_g3());
   lower->set_vd(upper->get_vd());
   lower->set_v(upper->get_v());
   lower->set_Lambdax(upper->get_Lambdax());
   lower->set_mu2(upper->get_mu2());

   const auto& Yv = upper->get_Yv();
   const auto& Mv = upper->get_Mv();
   lower->set_WOp(Yv.transpose() * Mv.inverse() * Yv);

   lower->set_scale(upper->get_scale());
}

double cSMHdCKM_cSMHdCKMRHN_matching::get_scale() const
{
    return scale;
}

double cSMHdCKM_cSMHdCKMRHN_matching::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

void cSMHdCKM_cSMHdCKMRHN_matching::set_models(Model *lower_, Model *upper_)
{
    lower = cast_model<cSMHdCKM<Two_scale>*>(lower_);
    upper = cast_model<cSMHdCKMRHN<Two_scale>*>(upper_);
}

} // namespace flexiblesusy
