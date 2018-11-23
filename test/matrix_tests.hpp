#ifndef MATRIX_TESTS_H
#define MATRIX_TESTS_H

#include "logger.hpp"
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
                 typename Eigen::NumTraits<Scalar>::Real max_deviation
                 = std::numeric_limits<typename Eigen::NumTraits<Scalar>::Real>::epsilon())
{
   const int rows = m.rows();
   const int cols = m.cols();

   bool result = true;
   for (int j = 0; j < cols; ++j) {
      for (int i = 0; i < rows; ++i) {
         if (i == j) {
            continue;
         }
         if (!is_zero(Abs(m(i,j)), max_deviation)) {
            result = false;
            break;
         }
      }
   }

   if (!result) {
      WARNING("matrix is not diagonal!"
              "\n   Matrix = " << m <<
              "\n   Maximum allowed max_deviation from zero = " << max_deviation);
   }

   return result;
}

template <typename Scalar, int N>
bool is_orthogonal(const Eigen::Matrix<Scalar,N,N>& m,
                   typename Eigen::NumTraits<Scalar>::Real max_deviation
                   = std::numeric_limits<typename Eigen::NumTraits<Scalar>::Real>::epsilon())
{
   const Eigen::Matrix<Scalar,N,N> unity(m * m.transpose());

   const bool result = is_equal(unity, Eigen::Matrix<Scalar,N,N>::Identity(), max_deviation);

   if (!result) {
      WARNING("matrix is not orthogonal!"
              "\n   Matrix = " << m <<
              "\n   Matrix * Matrix^T = " << unity <<
              "\n   Maximum allowed max_deviation from unity = " << max_deviation);
   }

   return result;
}

template <typename Scalar, int N>
bool is_unitary(const Eigen::Matrix<Scalar,N,N>& m,
                typename Eigen::NumTraits<Scalar>::Real max_deviation
                = std::numeric_limits<typename Eigen::NumTraits<Scalar>::Real>::epsilon())
{
   const Eigen::Matrix<Scalar,N,N> unity(m * m.adjoint());

   const bool result = is_equal(unity, Eigen::Matrix<Scalar,N,N>::Identity(), max_deviation);

   if (!result) {
      WARNING("matrix is not unitary!"
              "\n   Matrix = " << m <<
              "\n   Matrix * Matrix^+ = " << unity <<
              "\n   Maximum allowed max_deviation from unity = " << max_deviation);
   }

   return result;
}

template <typename Scalar, int N>
bool is_symmetric(const Eigen::Matrix<Scalar,N,N>& m,
                  typename Eigen::NumTraits<Scalar>::Real max_deviation
                  = std::numeric_limits<typename Eigen::NumTraits<Scalar>::Real>::epsilon())
{
   const bool result = is_equal(m, m.transpose(), max_deviation);

   if (!result) {
      WARNING("matrix is not symmetric!"
              "\n   Matrix = " << m <<
              "\n   Matrix^T = " << m.transpose() <<
              "\n   Maximum allowed max_deviation from equality = " << max_deviation);
   }

   return result;
}

} // namespace flexiblesusy

#endif
