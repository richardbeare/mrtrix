/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 27/06/08.

    This file is part of MRtrix.

    MRtrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRtrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef __math_cholesky_h__
#define __math_cholesky_h__

#include <gsl/gsl_linalg.h>
#include "math/matrix.h"

namespace MR {
  namespace Math {
    namespace Cholesky {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
      // float definitions of GSL functions:
      int gsl_linalg_cholesky_decomp (gsl_matrix_float * A);
      int gsl_linalg_cholesky_solve (const gsl_matrix_float * cholesky, const gsl_vector_float * b, gsl_vector_float * x);
      int gsl_linalg_cholesky_svx (const gsl_matrix_float * cholesky, gsl_vector_float * x);
      int gsl_linalg_cholesky_invert (gsl_matrix_float * cholesky);
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

      /** @addtogroup linalg 
        @{ */

      /** @defgroup cholesky Cholesky decomposition
        @{ */


      //! %Cholesky decomposition of \a A into lower triangular matrix
      /** \note the contents of \a A will be overwritten with its %Cholesky decomposition */
      template <typename T> inline MatrixView<T>& decomp (MatrixView<T>& A) {
       	gsl_linalg_cholesky_decomp (&A); 
	return (A);
      }

      //! solve A*x = b given its %Cholesky decomposition \a D 
      template <typename T> inline VectorView<T>& solve (VectorView<T>& x, const MatrixView<T>& D, const VectorView<T>& b) {
        gsl_linalg_cholesky_solve (&D, &b, &x);
        return (x);
      }

      //! solve A*x = b given its %Cholesky decomposition \a D, in place
      /** \note b is passed in as \a x */
      template <typename T> inline VectorView<T>& solve (VectorView<T>& x, const MatrixView<T>& D) {
        gsl_linalg_cholesky_svx (&D, &x);
        return (x);
      }

      //! invert A given its %Cholesky decomposition \a D, in place.
      template <typename T> inline MatrixView<T>& inv_from_decomp (MatrixView<T>& D) {
	gsl_linalg_cholesky_invert (&D);
	return (D);
      }

      //! invert \a A using %Cholesky decomposition, in place.
      template <typename T> inline MatrixView<T>& inv (MatrixView<T>& A) { return (inv_from_decomp (decomp (A))); }

      /** @} */
      /** @} */

    }
  }
}

#endif







