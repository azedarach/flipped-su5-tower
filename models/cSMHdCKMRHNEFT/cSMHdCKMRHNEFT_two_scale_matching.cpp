#include "cSMHdCKMRHNEFT_two_scale_matching.hpp"
#include "cSMHdCKMRHNEFT_matching.hpp"
#include "cSMHdCKMRHN_two_scale_model.hpp"
#include "cSMHdCKM_two_scale_model.hpp"

#include "error.hpp"
#include "single_scale_constraint.hpp"

// @todo remove
#include <iostream>

namespace flexiblesusy {

#define CLASSNAME cSMHdCKMRHNEFT_matching_up<Two_scale>

CLASSNAME::cSMHdCKMRHNEFT_matching_up(
   cSMHdCKM<Two_scale>* low_,
   cSMHdCKMRHN<Two_scale>* high_,
   const Scale_getter& scale_getter_,
   int loop_order_)
   : model(high_)
   , eft(low_)
   , scale_getter(scale_getter_)
   , scale(0.)
{
   set_loop_order(loop_order_);
}

double CLASSNAME::get_scale() const
{
   if (scale != 0.)
      return scale;

   if (!scale_getter)
      throw SetupError("Scale getter is not set!");

   return scale_getter();
}

void CLASSNAME::set_models(Model* low, Model* high)
{
   eft = cast_model<cSMHdCKM<Two_scale>*>(low);
   model = cast_model<cSMHdCKMRHN<Two_scale>*>(high);
}

void CLASSNAME::match()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   if (model->get_thresholds() && loop_order)
      cSMHdCKMRHNEFT_matching::match_low_to_high_scale_model(*model, *eft, loop_order);
   else
      cSMHdCKMRHNEFT_matching::match_low_to_high_scale_model_tree_level(*model, *eft);
}

void CLASSNAME::match_tree_level()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   cSMHdCKMRHNEFT_matching::match_low_to_high_scale_model_tree_level(*model, *eft);
}

int CLASSNAME::get_loop_order() const
{
   return loop_order;
}

void CLASSNAME::set_loop_order(int loop_order_)
{
   loop_order = loop_order_;
}

void CLASSNAME::set_scale(const Scale_getter& scale_getter_)
{
   scale_getter = scale_getter_;
}

void CLASSNAME::set_scale(Scale_getter&& scale_getter_)
{
   scale_getter = std::move(scale_getter_);
}

void CLASSNAME::set_scale(double scale_)
{
   scale = scale_;
}

#ifdef CLASSNAME
#undef CLASSNAME
#endif

#define CLASSNAME cSMHdCKMRHNEFT_matching_down<Two_scale>

CLASSNAME::cSMHdCKMRHNEFT_matching_down(
   cSMHdCKM<Two_scale>* low_,
   cSMHdCKMRHN<Two_scale>* high_,
   const Scale_getter& scale_getter_,
   int loop_order_)
   : model(high_)
   , eft(low_)
   , scale_getter(scale_getter_)
   , scale(0.)
{
   set_loop_order(loop_order_);
}

double CLASSNAME::get_scale() const
{
   if (scale != 0.)
      return scale;

   if (!scale_getter)
      throw SetupError("Scale getter is not set!");

   return scale_getter();
}

void CLASSNAME::set_models(Model* high, Model* low)
{
   eft = cast_model<cSMHdCKM<Two_scale>*>(low);
   model = cast_model<cSMHdCKMRHN<Two_scale>*>(high);
}

void CLASSNAME::match()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   std::cout << "start of downwards matching:\n";
   eft->print(std::cout);
   model->print(std::cout);

   if (model->get_thresholds() && loop_order)
      cSMHdCKMRHNEFT_matching::match_high_to_low_scale_model(*eft, *model, loop_order);
   else
      cSMHdCKMRHNEFT_matching::match_high_to_low_scale_model_tree_level(*eft, *model);

   std::cout << "end of downwards matching:\n";
   eft->print(std::cout);
   model->print(std::cout);
}

void CLASSNAME::match_tree_level()
{
   if (!model || !eft)
      throw SetupError("Model pointer in matching class is NULL!");

   eft->run_to(get_scale());
   model->run_to(get_scale());

   cSMHdCKMRHNEFT_matching::match_high_to_low_scale_model_tree_level(*eft, *model);
}

int CLASSNAME::get_loop_order() const
{
   return loop_order;
}

void CLASSNAME::set_loop_order(int loop_order_)
{
   loop_order = loop_order_;
}

void CLASSNAME::set_scale(const Scale_getter& scale_getter_)
{
   scale_getter = scale_getter_;
}

void CLASSNAME::set_scale(Scale_getter&& scale_getter_)
{
   scale_getter = std::move(scale_getter_);
}

void CLASSNAME::set_scale(double scale_)
{
   scale = scale_;
}

} // namespace flexiblesusy
