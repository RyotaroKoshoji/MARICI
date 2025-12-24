#include <iostream>

#include <mpi.h>

#include "MpiPolicy.h"
#include "IException.h"

#include "Program.h"


int main(int argc, char* argv[])
{
	int mpiRank = 0;
	int mpiProcessing = 1;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &mpiProcessing);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

	try
	{
		System::Parallel::MpiPolicy::setMpiRank(static_cast<std::size_t>(mpiRank));
		System::Parallel::MpiPolicy::setMpiProcessing(static_cast<std::size_t>(mpiProcessing));

		System::Marici::Program program(argc, argv);
		program.execute();
	}

	catch (const System::ExceptionServices::IException& e)
	{
		std::cerr << "Exception: ";
		std::cerr << e.toString() << "\n";
		std::cerr << "Terminate this program.\n\n";
	}

	catch (const std::exception& e)
	{
		std::cerr << "Exception: ";
		std::cerr << e.what() << "\n";
		std::cerr << "Terminate this program.\n\n";
	}


	std::cout << "MPI " << mpiRank << ": Finalize" << std::endl;
	MPI_Finalize();

	return 0;
}
