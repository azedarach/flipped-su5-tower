#ifndef MATRIX_TESTS_H
#define MATRIX_TESTS_H

#include "numerics2.hpp"
#include "wrappers.hpp"

#include <Eigen/Core>

#include <limits>
#include <type_traits>

namespace flexiblesusy {

template <class DerivedA, class DerivedB>
bool is_equal(const Eigen::MatrixBase<DerivedA>& a,
              const Eigen::MatrixBase<DerivedB>& b,
              typename std::common_type<typename Eigen::NumTraits<typename Eigen::MatrixBase<DerivedA>::Scalar>::Real,
              typename Eigen::NumTraits<typename Eigen::MatrixBase<DerivedB>::Scalar>::Real>::type max_dev)
{
   return (a - b).cwiseAbs().maxCoeff() <= max_dev;
}

template <typename Scalar, int N>
bool is_diagonal(const Eigen::Matrix<Scalar,N,N>& m,
                 typename Eigen::NumTraits<Scalar>::Real tolerance
                 = std::numeric_limits<typename Eigen::NumTraits<Scalar>::Real>::epsilon())
{
   const int rows = m.rows();
   const int cols = m.cols();

   for (int j = 0; j < cols; ++j) {
      for (int i = 0; i < rows; ++i) {
         if (i == j) {
            continue;
         }
         if (!is_zero(Abs(m(i,j)), tolerance)) {
            return false;
         }
      }
   }

   return true;
}

template <typename Scalar, int N>
bool is_orthogonal(const Eigen::Matrix<Scalar,N,N>& m,
                   typename Eigen::NumTraits<Scalar>::Real max_deviation
                   = std::numeric_limits<typename Eigen::NumTraits<Scalar>::Real>::epsilon())
{
   const Eigen::Matrix<Scalar,N,N> unity(m * m.transpose());

   return is_equal(unity, Eigen::Matrix<Scalar,N,N>::Identity(), max_deviation);
}

template <typename Scalar, int N>
bool is_unitary(const Eigen::Matrix<Scalar,N,N>& m,
                typename Eigen::NumTraits<Scalar>::Real max_deviation
                = std::numeric_limits<typename Eigen::NumTraits<Scalar>::Real>::epsilon())
{
   const Eigen::Matrix<Scalar,N,N> unity(m * m.adjoint());

   return is_equal(unity, Eigen::Matrix<Scalar,N,N>::Identity(), max_deviation);
}

template <typename Scalar, int N>
bool is_symmetric(const Eigen::Matrix<Scalar,N,N>& m,
                  typename Eigen::NumTraits<Scalar>::Real max_deviation
                  = std::numeric_limits<typename Eigen::NumTraits<Scalar>::Real>::epsilon())
{
   return is_equal(m, m.transpose(), max_deviation);
}

} // namespace flexiblesusy

#endif
