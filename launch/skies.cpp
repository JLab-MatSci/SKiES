#include <vector>
#include <string>
#include <unordered_map>

#include <launch/timer.h>
#include <launch/parser.h>
#include <launch/commands.h>

#include <skies/quantities/eigenvals.h>
#include <skies/lattices/latt_protocol.h>

#include <iostream>

#ifdef __cplusplus
extern "C" {
#endif

void epiInit();
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

int main(int argc, char **argv)
try {
	skies::launch::Timer t;
	t.start();

	print_name();
#ifdef SKIES_TBB
	std::cout << "TBB version launched" << std::endl << std::endl;
#else
	std::cout << "No TBB version launched" << std::endl << std::endl;
#endif
	t.print_start("The program started at ");

	std::vector<std::string> args;
	std::unordered_map<std::string, std::string> opts;
	std::string cmd = "help";
	if (argc > 1) cmd = argv[1];

	if (cmd != "help" && cmd != "list")
	{
		std::cout << "======== Some standard output from QE EPW:" << std::endl;
		epiInit();
	}
	skies::launch::parse_opts(argc - 1, argv + 1, args, opts);
	if ((cmd == "help") && (args.size() == 1))
		skies::launch::help_for_cmd(skies::launch::str_to_CMD(args[0]));
	skies::launch::resolve_cmd(skies::launch::str_to_CMD(cmd), args, opts);

	t.stop();
	t.print_stop("SKiES finished successfully at ");
	std::cout << "Total work time: " << t.elapsed() << " ms" << std::endl;

} catch (std::runtime_error& err) {
	std::cout << "\nERROR: " << err.what() << std::endl;
} catch (std::invalid_argument& err) {
	std::cout << "\nERROR: one or more given parameters are incorrect, error in function " << err.what() << std::endl;
}