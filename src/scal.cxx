#include <scalapackpp/wrappers/pblas/scal.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

void psscal_( const scalapack_int* N, const float* ALPHA, float* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, const scalapack_int* INCX );
void pdscal_( const scalapack_int* N, const double* ALPHA, double* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, const scalapack_int* INCX );
void pcscal_( const scalapack_int* N, const scomplex* ALPHA, scomplex* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, const scalapack_int* INCX );
void pzscal_( const scalapack_int* N, const dcomplex* ALPHA, dcomplex* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, const scalapack_int* INCX );

}


namespace scalapackpp::wrappers {

template <>
void pscal( scalapack_int N, float ALPHA, float * X, scalapack_int IX, scalapack_int JX, 
            const scalapack_desc& DESCX, scalapack_int INCX ) {

  psscal_( &N, &ALPHA, X, &IX, &JX, DESCX.data(), &INCX );

}
template <>
void pscal( scalapack_int N, double ALPHA, double * X, scalapack_int IX, scalapack_int JX, 
            const scalapack_desc& DESCX, scalapack_int INCX ) {

  pdscal_( &N, &ALPHA, X, &IX, &JX, DESCX.data(), &INCX );

}
template <>
void pscal( scalapack_int N, scomplex ALPHA, scomplex * X, scalapack_int IX, scalapack_int JX, 
            const scalapack_desc& DESCX, scalapack_int INCX ) {

  pcscal_( &N, &ALPHA, X, &IX, &JX, DESCX.data(), &INCX );

}
template <>
void pscal( scalapack_int N, dcomplex ALPHA, dcomplex * X, scalapack_int IX, scalapack_int JX, 
            const scalapack_desc& DESCX, scalapack_int INCX ) {

  pzscal_( &N, &ALPHA, X, &IX, &JX, DESCX.data(), &INCX );

}

}
