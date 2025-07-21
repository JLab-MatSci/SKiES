/*-----------------------------------------------------------------------
    * SKiES - Solver of Kinetic Equation for Solids
    * 
    * (C) 2025 Galtsov Ilya, Fokin Vladimir, Minakov Dmitry, Levashov Pavel (JIHT RAS)
    *
    * SKiES may only be utilized for non-profit research.
    * Citing appropriate sources is required when using SKiES.
    * 
    * Distribution of this file is permitted by the GNU General Public License.
    * Examine the `LICENSE' file located in the current distribution's root directory.
------------------------------------------------------------------------- */
#include <vector>
#include <string>
#include <unordered_map>

#include <launch/timer.h>
#include <launch/parser.h>
#include <launch/commands.h>

#include <skies/quantities/eigenvals.h>
#include <skies/lattices/latt_protocol.h>

#include <skies/utils/mpi_wrapper.h>

#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif

void epiInit(int*, int*);
void epiFinalize();

MPI_Fint MPICommC2F(MPI_Comm *comm) {
	MPI_Fint ret = MPI_Comm_c2f(*comm);
	return ret;
}

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
	std::cout << "SKiES software v.1.0.0\n(C) 2025 Galtsov Ilya, Fokin Vladimir, Minakov Dmitry, Levashov Pavel (JIHT RAS)\n\n";
}

int main(int argc, char **argv)
try 
{
	using namespace skies::utils;
	mpi::mpi_handler handler;

	skies::launch::Timer t;
	if (mpi::is_root()) {
		std::cout << "rank = " << mpi::rank() << std::endl;
		t.start();

		print_name();
#ifdef SKIES_TBB
		std::cout << "TBB version launched" << std::endl << std::endl;
#else
		std::cout << "No TBB version launched" << std::endl << std::endl;
#endif

		t.print_start("The program started at ");
	}

	std::vector<std::string> args;
	std::unordered_map<std::string, std::string> opts;
	std::string cmd = "help";
	if (argc > 1) cmd = argv[1];

	if (cmd != "help" && cmd != "list")
	{
		if (mpi::is_root())
			std::cout << "======== Some standard output from QE EPW:" << std::endl;
		int rank = mpi::rank();
		MPI_Fint f_comm = 0;
		epiInit(&f_comm, &rank);
	}
	skies::launch::parse_opts(argc - 1, argv + 1, args, opts);
	if ((mpi::is_root()) && (cmd == "help") && (args.size() == 1))
		skies::launch::help_for_cmd(skies::launch::str_to_CMD(args[0]));
	skies::launch::resolve_cmd(skies::launch::str_to_CMD(cmd), args, opts);

	if (mpi::is_root()) {
		t.stop();
		t.print_stop("SKiES finished successfully at ");
		std::cout << "Total work time: " << t.elapsed() << " s" << std::endl;
	}
} catch (std::runtime_error& err) {
	std::cout << "\nERROR: " << err.what() << std::endl;
} catch (std::invalid_argument& err) {
	std::cout << "\nERROR: one or more given parameters are incorrect, error in function " << err.what() << std::endl;
}
