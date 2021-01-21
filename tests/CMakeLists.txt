include( HandleCatch2 )

# BLAS++/LAPACK++
include( FetchContent )

if( Threads::Threads IN_LIST ScaLAPACK_LIBRARIES )
  find_package( Threads )
endif()

if( OpenMP::OpenMP_C IN_LIST ScaLAPACK_LIBRARIES )
  find_package( OpenMP )
endif()

if( MPI::MPI_C IN_LIST ScaLAPACK_LIBRARIES )
  find_package( MPI )
endif()

set( BLAS_LIBRARIES   "${ScaLAPACK_LIBRARIES}" )
set( LAPACK_LIBRARIES "${ScaLAPACK_LIBRARIES}" )

message( STATUS "BLAS_LIBRARIES = ${BLAS_LIBRARIES}" )
message( STATUS "LAPACK_LIBRARIES = ${LAPACK_LIBRARIES}" )

set( use_openmp FALSE CACHE BOOL "Disable OpenMP in ICL Wrappers" )

FetchContent_Declare( blaspp
  GIT_REPOSITORY https://bitbucket.org/icl/blaspp.git 
)
FetchContent_Declare( lapackpp
  GIT_REPOSITORY https://bitbucket.org/icl/lapackpp.git 
)

FetchContent_MakeAvailable( blaspp )
FetchContent_MakeAvailable( lapackpp )
target_compile_definitions( lapackpp PUBLIC LAPACK_COMPLEX_CPP )


add_library( ut_framework ut.cxx )
target_link_libraries( ut_framework PUBLIC scalapackpp scalapackpp::catch2 lapackpp )

add_executable( test_scalapackpp 
  scatter_gather.cxx 
  gemm.cxx 
  trmm.cxx 
  trsm.cxx 
  potrf.cxx 
  sevp.cxx 
  index.cxx 
  geadd.cxx 
  gesvd.cxx
  transpose.cxx 
  gevp.cxx 
  potrs.cxx 
  posv.cxx 
  gesv.cxx 
  getrf_s.cxx 
  trtrs.cxx 
  potri.cxx 
  getri.cxx
  trtri.cxx
  norm.cxx
)
target_link_libraries( test_scalapackpp PUBLIC ut_framework )

add_test( NAME SCALAPACKPP_TEST 
          COMMAND
            ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_scalapackpp> ${MPIEXEC_POSTFLAGS}
)
