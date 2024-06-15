#include <vector>
#include <string>
#include <unordered_map>

#include <launch/timer.h>
#include <launch/parser.h>
#include <launch/commands.h>

#include <skies/utils/mpi_wrapper.h>
#include <skies/quantities/eigenvals.h>
#include <skies/lattices/latt_protocol.h>

#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif

void epiInit(int* rank);
void epiFinalize();

#ifdef __cplusplus
}
#endif

void print_name()
{
	std::cout <<" _____ _   ___ _____ _____" << std::endl;
	std::cout <<"/  ___| | / (_)  ___/  ___|" << std::endl;
	std::cout << "\\ `--.| |/ / _| |__ \\ `--." << std::endl;
	std::cout << " `--. \\    \\| |  __| `--. \\" << std::endl;
	std::cout << "/\\__/ / |\\  \\ | |___/\\__/ /" << std::endl;
	std::cout << "\\____/\\_| \\_/_\\____/\\____/" << std::endl;
	std::cout << std::endl;
	std::cout << "SKiES software v.1.0.0\n(C) 2024 Galtsov Ilya, Minakov Dmitry, Fokin Vladimir, Levashov Pavel (JIHT RAS)\n\n";
}

void print_numproc(int nproc)
{
#ifdef SKIES_MPI
	std::cout << "Parallel MPI version. Launched on " << nproc  << " processors" << std::endl;
#else
	std::cout << "Serial version launched" << std::endl << std::endl;
#endif
}

int main(int argc, char **argv)
try {
#ifdef SKIES_MPI
	skies::mpi::init(&argc, &argv);
#endif
	int rank{ 0 };
	int nproc{ 1 };
#ifdef SKIES_MPI
	rank = skies::mpi::rank();
	nproc = skies::mpi::size();
#endif
	skies::launch::Timer t;
	t.start();

    if (!rank)
	{
		print_name();
		print_numproc(nproc);
		t.print_start("The program started at ");
	}

	std::vector<std::string> args;
	std::unordered_map<std::string, std::string> opts;
	std::string cmd = "help";
	if (argc > 1) cmd = argv[1];

	if (cmd != "help" && cmd != "list")
	{
		if (!rank) std::cout << "======== Some standard output from QE EPW:" << std::endl;
		epiInit(&rank);
	}
	skies::launch::parse_opts(argc - 1, argv + 1, args, opts);
	skies::launch::exec_cmd(cmd, args, opts);

	t.stop();
	if (!rank)
	{
		t.print_stop("SKiES finished successfully at ");
		std::cout << "Total work time: " << t.elapsed() << " ms" << std::endl;
	}

#ifdef SKIES_MPI
	if (cmd != "help" && cmd != "list")
		epiFinalize();
	skies::mpi::finalize();
#endif
} catch (std::runtime_error& err) {
	std::cout << "\nERROR: " << err.what() << std::endl;
} catch (std::invalid_argument& err) {
	std::cout << "\nERROR: one or more given parameters are incorrect, error in function " << err.what() << std::endl;
}