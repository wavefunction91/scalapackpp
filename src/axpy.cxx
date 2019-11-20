#include <scalapackpp/wrappers/pblas/axpy.hpp>

using scalapackpp::scalapack_int;
using scalapackpp::dcomplex;
using scalapackpp::scomplex;

// Prototypes
extern "C" {

void psaxpy_( const scalapack_int* N, const float* ALPHA, const float* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, float* Y, const scalapack_int* IY, 
              const scalapack_int* JY, const scalapack_int* DESCY );
void pdaxpy_( const scalapack_int* N, const double* ALPHA, const double* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, double* Y, const scalapack_int* IY, 
              const scalapack_int* JY, const scalapack_int* DESCY );
void pcaxpy_( const scalapack_int* N, const scomplex* ALPHA, const scomplex* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, scomplex* Y, const scalapack_int* IY, 
              const scalapack_int* JY, const scalapack_int* DESCY );
void pzaxpy_( const scalapack_int* N, const dcomplex* ALPHA, const dcomplex* X, 
              const scalapack_int* IX, const scalapack_int* JX, 
              const scalapack_int* DESCX, dcomplex* Y, const scalapack_int* IY, 
              const scalapack_int* JY, const scalapack_int* DESCY );

}


namespace scalapackpp::wrappers {

template <>
void paxpy( scalapack_int N, float ALPHA,
    const float * X, scalapack_int IX, scalapack_int JX, const scalapack_desc& DESCX,
          float * Y, scalapack_int IY, scalapack_int JY, const scalapack_desc& DESCY ){

  psaxpy_( &N, &ALPHA, X, &IX, &JX, DESCX.data(), Y, &IY, &JY, DESCY.data() );

}
template <>
void paxpy( scalapack_int N, double ALPHA,
    const double * X, scalapack_int IX, scalapack_int JX, const scalapack_desc& DESCX,
          double * Y, scalapack_int IY, scalapack_int JY, const scalapack_desc& DESCY ){

  pdaxpy_( &N, &ALPHA, X, &IX, &JX, DESCX.data(), Y, &IY, &JY, DESCY.data() );

}
template <>
void paxpy( scalapack_int N, scomplex ALPHA,
    const scomplex * X, scalapack_int IX, scalapack_int JX, const scalapack_desc& DESCX,
          scomplex * Y, scalapack_int IY, scalapack_int JY, const scalapack_desc& DESCY ){

  pcaxpy_( &N, &ALPHA, X, &IX, &JX, DESCX.data(), Y, &IY, &JY, DESCY.data() );

}
template <>
void paxpy( scalapack_int N, dcomplex ALPHA,
    const dcomplex * X, scalapack_int IX, scalapack_int JX, const scalapack_desc& DESCX,
          dcomplex * Y, scalapack_int IY, scalapack_int JY, const scalapack_desc& DESCY ){

  pzaxpy_( &N, &ALPHA, X, &IX, &JX, DESCX.data(), Y, &IY, &JY, DESCY.data() );

}

}
