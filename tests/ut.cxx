#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#include <mpi.h>

int main(int argc, char* argv[])
{
    MPI_Init(&argc,&argv);

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    int res = Catch::Session().run(argc, argv);
    MPI_Finalize();

    return res;
}
