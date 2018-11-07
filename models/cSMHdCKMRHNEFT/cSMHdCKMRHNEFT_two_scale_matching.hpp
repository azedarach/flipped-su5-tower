#ifndef cSMHdCKMRHNEFT_TWO_SCALE_MATCHING_H
#define cSMHdCKMRHNEFT_TWO_SCALE_MATCHING_H

#include "single_scale_matching.hpp"

#include <functional>

namespace flexiblesusy {

class Model;
class Two_scale;

template <class T> class cSMHdCKM;
template <class T> class cSMHdCKMRHN;
template <class T> class cSMHdCKMRHNEFT_matching_up;
template <class T> class cSMHdCKMRHNEFT_matching_down;

template<>
class cSMHdCKMRHNEFT_matching_up<Two_scale> : public Single_scale_matching {
public:
   using Scale_getter = std::function<double()>;

   cSMHdCKMRHNEFT_matching_up() = default;

   cSMHdCKMRHNEFT_matching_up(cSMHdCKM<Two_scale>*,
                              cSMHdCKMRHN<Two_scale>*,
                              const Scale_getter&,
                              int);

   virtual ~cSMHdCKMRHNEFT_matching_up() = default;

   virtual double get_scale() const override;
   virtual void set_models(Model*, Model*) override;
   virtual void match() override;

   int get_loop_order() const;
   void set_loop_order(int);
   void set_scale(double);
   void set_scale(const Scale_getter&);
   void set_scale(Scale_getter&&);
   void match_tree_level();

private:
   cSMHdCKMRHN<Two_scale>* model{nullptr};
   cSMHdCKM<Two_scale>* eft{nullptr};
   Scale_getter scale_getter{};
   double scale{0.};
   int loop_order{0};
};

template<>
class cSMHdCKMRHNEFT_matching_down<Two_scale> : public Single_scale_matching {
public:
   using Scale_getter = std::function<double()>;

   cSMHdCKMRHNEFT_matching_down() = default;

   cSMHdCKMRHNEFT_matching_down(cSMHdCKM<Two_scale>*,
                                cSMHdCKMRHN<Two_scale>*,
                                const Scale_getter&,
                                int);

   virtual ~cSMHdCKMRHNEFT_matching_down() = default;

   virtual double get_scale() const override;
   virtual void set_models(Model*, Model*) override;
   virtual void match() override;

   int get_loop_order() const;
   void set_loop_order(int);
   void set_scale(double);
   void set_scale(const Scale_getter&);
   void set_scale(Scale_getter&&);
   void match_tree_level();

private:
   cSMHdCKMRHN<Two_scale>* model{nullptr};
   cSMHdCKM<Two_scale>* eft{nullptr};
   Scale_getter scale_getter{};
   double scale{0.};
   int loop_order{0};
};

} // namespace flexiblesusy

#endif
