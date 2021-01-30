#include <scalapackpp/block_cyclic.hpp>

int main() {
  MPI_Init( NULL, NULL );
  {
    auto grid = std::make_shared<const blacspp::Grid>(blacspp::Grid::square_grid( MPI_COMM_WORLD ));
    scalapackpp::BlockCyclicDist2D dist( grid, 2, 2 );
  }
  MPI_Finalize();
}
