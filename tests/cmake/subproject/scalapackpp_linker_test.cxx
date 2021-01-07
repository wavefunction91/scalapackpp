#include <scalapackpp/block_cyclic.hpp>

int main() {
  MPI_Init( NULL, NULL );
  {
    blacspp::Grid grid = blacspp::Grid::square_grid( MPI_COMM_WORLD );
    scalapackpp::BlockCyclicDist2D dist( grid, 2, 2 );
  }
  MPI_Finalize();
}
