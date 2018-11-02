#ifndef cSMHdCKM_cSMHdCKMRHN_TWO_SCALE_MATCHING_H
#define cSMHdCKM_cSMHdCKMRHN_TWO_SCALE_MATCHING_H

#include "cSMHdCKM_cSMHdCKMRHN_matching.hpp"
#include "single_scale_matching.hpp"

namespace flexiblesusy {

template <class T>
class cSMHdCKM;

template <class T>
class cSMHdCKMRHN;

class Two_scale;

class cSMHdCKM_cSMHdCKMRHN_matching {
public:

   void match_low_to_high_scale_model();
   void match_high_to_low_scale_model();
   double get_scale() const;
   void set_models(Model* lower, Model* upper);
   double get_initial_scale_guess() const;
   void set_scale(double) {}
   void reset();

private:
   double scale{0.};
   double initial_scale_guess{0.};
   cSMHdCKM<Two_scale>* lower{nullptr};
   cSMHdCKMRHN<Two_scale>* upper{nullptr};
};

template<>
class cSMHdCKM_cSMHdCKMRHN_matching_up<Two_scale>
   : public Single_scale_matching {
public:
   cSMHdCKM_cSMHdCKMRHN_matching_up();
   void match();
   double get_scale() const;
   void set_models(Model *lower, Model *upper);
   double get_initial_scale_guess() const;
   void set_scale(double);
   void reset();

private:
   cSMHdCKM_cSMHdCKMRHN_matching matching{};
};

template<>
class cSMHdCKM_cSMHdCKMRHN_matching_down<Two_scale>
   : public Single_scale_matching {
public:
   cSMHdCKM_cSMHdCKMRHN_matching_down();
   void match();
   double get_scale() const;
   void set_models(Model *lower, Model *upper);
   double get_initial_scale_guess() const;
   void set_scale(double);
   void reset();

private:
   cSMHdCKM_cSMHdCKMRHN_matching matching{};
};

} // namespace flexiblesusy

#endif
