#include <iostream>
#include <stdexcept>

#include <skies/common/alg.h>
#include <skies/common/units.h>
#ifdef SKIES_MPI
#include <skies/utils/mpi_wrapper.h>
#endif
#include <skies/common/ndimarrays.h>

#include <skies/quantities/dos.h>
#include <skies/quantities/gmatrix.h>
#include <skies/quantities/eigenvals.h>
#include <skies/quantities/highSymPath.h>

#include <skies/sampling/tetrahedra.h>

#include <skies/lattices/kp_protocol.h>
#include <skies/lattices/latt_protocol.h>

#include <skies/spectral/spec_func.h>

#include <skies/transport/elastic.h>
#include <skies/transport/inelastic.h>

#include <launch/timer.h>
#include <launch/commands.h>

namespace skies { namespace launch {

using namespace arrays;
using namespace spectral;
using namespace quantities;

namespace {
const std::string HELPMSG = \
"SKiES help message:\n"
"skies list                 the list of all supported commands\n"
"skies help [command]       the description of the chosen command\n"
"skies [command] [options]  launch the chosen command with desired options\n";
}

template <typename T>
std::vector<T> parse_vector_of_numbers(const std::string& opts, const std::string& name, size_t argc)
{
	assert(argc > 0);

	std::string numbers = "=[";
	for (size_t i = 0; i < argc - 1; ++i) numbers += "num,";
	numbers += "num]\n";
	const std::string error = "Incorrect vector option: use " + name + numbers;

	if (opts[0] != '[' || opts[opts.length() - 1] != ']')
		throw std::runtime_error(error);

	auto inner = opts.substr(1, opts.length() - 2);
	auto splitted_line = custom_split(inner, ',');
	if (splitted_line.size() != argc)
		throw std::runtime_error(error);

	std::vector<T> args;
	if (std::is_same_v<int, std::remove_reference_t<T>>)
	{
		for (size_t i = 0; i < argc; ++i)
			args.push_back(std::stoi(splitted_line[i].data()));
	}
	if (std::is_same_v<size_t, std::remove_reference_t<T>>)
	{
		for (size_t i = 0; i < argc; ++i)
			args.push_back(static_cast<size_t>(std::stoi(splitted_line[i].data())));
	}
	else if (std::is_same_v<double, std::remove_reference_t<T>>)
	{
		for (size_t i = 0; i < argc; ++i)
			args.push_back(std::stod(splitted_line[i].data()));
	}

	return args;
}

void exec_cmd(const std::string& cmd,
	         std::vector<std::string>& args,
	         std::unordered_map<std::string, std::string>& opts)
{
    int rank{ 0 };
#ifdef SKIES_MPI
	rank = mpi::rank();
#endif
	if (!resolve_cmd(cmd, args, opts) && !rank)
    {
        std::string wrong_cmd = "Sorry, command " + cmd + " is not supported. Please, use the help command";
		throw std::runtime_error(wrong_cmd);
    }
}

bool resolve_cmd(const std::string& cmd,
                 std::vector<std::string>& args,
                 std::unordered_map<std::string, std::string>& opts)
{
	bool cmd_found = false;
	int rank = 0;
#ifdef SKIES_MPI
	rank = mpi::rank();
#endif
	if (cmd == "help" || cmd == "list")
    {
		if(cmd == "list" || !args.size())
            cmd_found = true;
		if (!rank)
        {
			if (cmd == "help" && !args.size())
				std::cout << HELPMSG;
			if (cmd == "list")
				std::cout << "SKiES: supported commands list:\n";
		}
	}

	CMD("dos",
		"calculates electron density of states using Wannier interpolation (EPW)",
		"skies dos [options]:\n"
		" Mandatory:\n"
		"  --grid=<[nk1,nk2,nk3]>: Monchorst-Pack grid of k-points (nk1 x nk2 x nk3) in the 1st BZ\n"
		" Non-mandatory:\n"
		"  Options can be given in any order. Options include:\n"
		"  --eF=<num>: Fermi level (in eV). If given energy axis zero will be at given Fermi level\n"
        "  --range=<[low,high]>: electron energy range to bracket, from low to high (in eV). Default is [-10.0,10.0]\n"
        "  --sampling=<str>: type of smearing for delta-functions approximation. Supported values are: 'gs' (gaussians), 'fd' (fermi-dirac derivative). Default value is 'fd'\n"
		"  --smearing=<sigma>: value of smearing for chosen sampling in eV. Default value is 0.03\n"
        "  --bins=<num>: amount of bins. Default is 300\n"
		"  --tetra: if given, use tetrahedron method for BZ sampling\n"
	) {
		EigenValue::eF = opts["eF"] == "" ? 0 : std::stod(opts["eF"]);

		if (opts["grid"] == "")
			throw std::runtime_error("k-point grid is not specified");
        
		auto grid = parse_vector_of_numbers<size_t>(opts["grid"], "grid", 3);

        double low_energy{ -10.0 };
        double high_energy{ 10.0 };
		if (opts["range"] != "")
        {
			auto epsilons = parse_vector_of_numbers<double>(opts["range"], "range", 2);
			low_energy  = epsilons[0];
			high_energy = epsilons[1]; 
        }

        std::string sampling = "fd";
        if (opts["sampling"] != "")
		    sampling = opts["sampling"];

        double smearing = 0.03;
        if (opts["smearing"] != "")
		    smearing = stod(opts["smearing"]);

        size_t bins = 300;
        if (opts["bins"] != "")
		    bins = static_cast<size_t>(stoi(opts["bins"]));

		bool is_tetra = opts["tetra"] == "true" ? true : false;

		launch::Timer t;
		t.start("========= Evaluating electronic DOS using Wannier interpolation (EPW)...");
        
		KPprotocol kprot(grid[0], grid[1], grid[2]);
		tetrahedra::TetraHandler::set_kprot(kprot);
		array2D weights(kprot.nkpt(), array1D(EigenValue::nbands, 1));
		auto range = arrays::create_range(low_energy, high_energy, bins);

		if (is_tetra)
		{
			tetrahedra::evaluate_dos(range);
		}
		else
		{
        	quantities::evaluate_dos<EigenValueDrawable>(kprot.grid(), range, smearing, bzsampling::switch_sampling(sampling), weights);
		}

		t.stop("========= Electron DOS is evaluated. Results are written to EigenValueDOS.dat");
		t.print_elapsed("\t  DOS evaluation time: ");
	} CMD_END;

	CMD("phdos",
		"calculates phonon density of states using Wannier interpolation",
		"skies phdos [options]:\n"
		" Mandatory:\n"
		"  --grid=<[nk1,nk2,nk3]>: Monchorst-Pack grid of k-points (nk1 x nk2 x nk3) in the 1st BZ\n"
		" Non-mandatory:\n"
		"  Options can be given in any order. Options include:\n"
        "  --range=<[low,high]>: electron energy range to bracket, from low to high in meV. Default is [0.0, 50.0]\n"
        "  --sampling=<str>: type of smearing for delta-functions approximation. Supported values are: 'gs' (gaussians), 'fd' (fermi-dirac derivative). Default value is 'fd'\n"
		"  --smearing=<sigma>: value of smearing for chosen sampling in meV. Default value is 0.5\n"
        "  --bins=<num>: amount of bins. Default is 300\n"
		"  --tetra: if given, use tetrahedron method for BZ sampling\n"
	) {
		if (opts["grid"] == "")
			throw std::runtime_error("k-point grid is not specified");
        
		auto grid = parse_vector_of_numbers<size_t>(opts["grid"], "grid", 3);

        double low_energy{ 0.0 / 1000.0 };
        double high_energy{ 50.0 / 1000.0 };
		if (opts["range"] != "")
        {
			auto epsilons = parse_vector_of_numbers<double>(opts["range"], "range", 2);
			low_energy  = epsilons[0] / 1000.0;
			high_energy = epsilons[1] / 1000.0; 
        }

        std::string sampling = "fd";
        if (opts["sampling"] != "")
		    sampling = opts["sampling"];

        double smearing = 0.0005;
        if (opts["smearing"] != "")
		    smearing = stod(opts["smearing"]) / 1000.0;

        size_t bins = 300;
        if (opts["bins"] != "")
		    bins = static_cast<size_t>(stoi(opts["smearing"]));

		bool is_tetra = opts["tetra"] == "true" ? true : false;

		launch::Timer t;
		t.start("========= Evaluating phonon DOS using Wannier interpolation (EPW)...");
        KPprotocol kprot(grid[0], grid[1], grid[2]);
		tetrahedra::TetraHandler::set_kprot(kprot);
        array2D weights(kprot.nkpt(), array1D(EigenFrequency::nmodes, 1));
		auto range = arrays::create_range(low_energy, high_energy, bins);

		if (is_tetra)
		{
			tetrahedra::evaluate_phdos(range);
		}
		else
		{
        	quantities::evaluate_dos<EigenFrequencyDrawable>(kprot.grid(), range, smearing, bzsampling::switch_sampling(sampling), weights);
		}

		t.stop("========= Phonon DOS is evaluated. Results are written to EigenFrequencyDOS.dat");
		t.print_elapsed("\t  PhDOS evaluation time: ");
	} CMD_END;

	CMD("trdos",
		"calculates transport density of states",
		"skies trdos [options]:\n"
		" Mandatory:\n"
		"  --grid=<[nk1,nk2,nk3]>: Monchorst-Pack grid of k-points (nk1 x nk2 x nk3) in the 1st BZ\n"
		" Non-mandatory:\n"
		"  Options can be given in any order. Options include:\n"
		"  --eF=<num>: Fermi level (in eV). If given energy axis zero will be at given Fermi level\n"
        "  --range=<[low,high]>: electron energy range to bracket, from low to high (in eV, relative to eF). Default is [-2.0,2.0]\n"
        "  --sampling=<str>: type of smearing for delta-functions approximation. Supported values are: 'gs' (gaussians), 'fd' (fermi-dirac derivative). Default value is 'fd'\n"
		"  --smearing=<sigma>: value of smearing for chosen sampling in eV. Default value is 0.03\n"
        "  --bins=<num>: amount of bins. Default is 100\n"
		"  --smeared: if given calculates smeared transport DOS. Requires file 'VelocitiesDOS.dat' to be precalculated\n"
		"  --Te=<num>: used for smeared DOS calculation, defines smearing value in eV. Default is 0.258\n"
		// "  --cart=<x, y or z>: cartesian index of velocities. Default is x\n"
		"  --tetra: if given, use tetrahedron method for BZ sampling\n"
	) {
		EigenValue::eF = opts["eF"] == "" ? 0 : std::stod(opts["eF"]);
		bool is_smeared = opts["smeared"] == "true" ? true : false;
		double sigma = opts["Te"] == "" ? 0.258 : std::stod(opts["Te"]);
		if (is_smeared)
		{
			launch::Timer t;
			t.start("========= Evaluating smeared transport DOS...");
			std::ifstream ifs("VelocitiesDOS.dat");
			if (ifs.fail())
				throw std::runtime_error("Smeared transport DOS requires precalculated transport DOS in VelocitiesDOS.dat");
			std::string line;
			array1D eigenvals, transDOSes;
			getline(ifs, line);
			while (ifs.good())
			{
				getline(ifs, line);
				if (!line.empty()) {
					auto splitted_line = custom_split(line, ' ');
					eigenvals.push_back(std::stod(splitted_line[0].data()));
					transDOSes.push_back(std::stod(splitted_line[1].data()));
				}
			}
			ifs.close();
			assert(eigenvals.size() == transDOSes.size());

			array1D transDOSes_smeared(transDOSes.size(), 0.0);
			std::transform(eigenvals.begin(), eigenvals.end(), transDOSes_smeared.begin(),
						[transDOSes, eigenvals, sigma] (double e) { return bzsampling::smear_with_fd(transDOSes, eigenvals, e, sigma); });
			std::ofstream os("VelocitiesSmearedDOS.dat");
			os << std::right;
			os << std::setw(13) << "# Energy [eV]";
			os << std::setw(49) << "    Smeared transport DOS [13.605685 * Ry bohr^2]";
			os << std::endl;
			for (size_t i = 0; i < eigenvals.size(); ++i)
				os << std::setprecision(6) << std::setw(13) << eigenvals[i] << std::setw(49) << transDOSes_smeared[i] << std::endl;
			os.close();
			t.stop("========= Smeared transport DOS is evaluated. Results are written to VelocitiesSmearedDOS.dat");
			t.print_elapsed("\t Transport DOS evaluation time: ");
			return cmd_found;
		}

		if (opts["grid"] == "")
			throw std::runtime_error("k-point grid is not specified");
        
		auto grid = parse_vector_of_numbers<size_t>(opts["grid"], "grid", 3);

        double low_energy{ -10.0 };
        double high_energy{ 10.0 };
		if (opts["range"] != "")
        {
			auto epsilons = parse_vector_of_numbers<double>(opts["range"], "range", 2);
			low_energy  = epsilons[0];
			high_energy = epsilons[1]; 
        }

        std::string sampling = "fd";
        if (opts["sampling"] != "")
		    sampling = opts["sampling"];

        double smearing = 0.03;
        if (opts["smearing"] != "")
		    smearing = stod(opts["smearing"]);

        size_t bins = 100;
        if (opts["bins"] != "")
		    bins = static_cast<size_t>(stoi(opts["bins"]));

		// char cart = opts["cart"] == "" ? 'x' : *opts["cart"].c_str();

		bool is_tetra = opts["tetra"] == "true" ? true : false;

		launch::Timer t;
		t.start("========= Evaluating transport DOS using Wannier interpolation (EPW)...");
		KPprotocol kprot(grid[0], grid[1], grid[2]);
		tetrahedra::TetraHandler::set_kprot(kprot);
		auto range = arrays::create_range(low_energy, high_energy, bins);

		if (is_tetra)
		{
			tetrahedra::evaluate_trdos(range);
		}
		else
		{
			quantities::evaluate_trdos(kprot.grid(), range, smearing, bzsampling::switch_sampling(sampling));
		}

		t.stop("========= Transport DOS is evaluated. Results are written to VelocitiesDOS.dat");
		t.print_elapsed("\t Transport DOS evaluation time: ");
	} CMD_END;

	CMD("bands",
		"calculates electron band structure along given k-path",
		"skies bands [options]:\n"
		"  Options can be given in any order. Options include:\n"
		"  --bins=<num>: number of bins between each pair of k-points in given path\n"
		"  --infile=<str>: the name of file with enumerated k-points.\n\n"
		"                  The format must be:\n"
		"                  <k-point-1 label> <num> <num> <num>\n"
		"                  <k-point-2 label> <num> <num> <num>\n"
		"                  ...................................\n"
		"                  <k-point-N label> <num> <num> <num>\n"
		"                  Here <num> after labels are crystal coordinates.\n"
		"                  Default name is 'KPATH'\n"
		"   --eF=<num>: if given band structure energies will be centered at given Fermi level\n"
	) {
		EigenValue::eF = opts["eF"] == "" ? 0 : std::stod(opts["eF"]);

		int bins = opts["bins"] == "" ? 30 : std::stoi(opts["bins"]);
		std::string	infile = opts["infile"] == "" ? "KPATH" : opts["infile"];

		std::ifstream ifs(infile);
		if (ifs.fail())
			throw std::runtime_error("The input file does not exist.\n");\

		interpol::kpLabelCoords kpath;
		std::string line;
		launch::Timer t;
		t.start("========= Evaluating band structure...");
		while (ifs.good())
		{
			getline(ifs, line);
			if (!line.empty())
			{
				auto line_splitted = custom_split(line, ' ');
				const std::string label = line_splitted[0].data();
				double kx = std::stod(line_splitted[1].data());
				double ky = std::stod(line_splitted[2].data());
				double kz = std::stod(line_splitted[3].data());
				
				array1D k{ kx, ky, kz };
				kpath.push_back(std::make_pair(label, k));
			}
		}
		interpol::evaluate_along_kp_path<EigenValueDrawable>(kpath, bins);
		t.stop("========= Band structure is calculated. The results are written to EigenValue.dat.\n"
			   "          See also KLABELS for correspondance between k-point labels and k-path distance.\n");
		t.print_elapsed("\t  Band structure evaluation time: ");
		ifs.close();
	} CMD_END;

	CMD("phonons",
		"calculates phonon dispersion curve along given k-path",
		"skies phonons [options]:\n"
		"  Options can be given in any order. Options include:\n"
		"  --bins=<num>: number of bins between each pair of k-points in given path\n"
		"  --infile=<str>: the name of file with enumerated k-points.\n\n"
		"                  The format must be:\n"
		"                  <k-point-1 label> <num> <num> <num>\n"
		"                  <k-point-2 label> <num> <num> <num>\n"
		"                  ...................................\n"
		"                  <k-point-N label> <num> <num> <num>\n"
		"                  Here <num> after labels are crystal coordinates.\n"
		"                  Default name is 'KPATH'\n"
	) {
		int bins = opts["bins"] == "" ? 30 : std::stoi(opts["bins"]);
		std::string	infile = opts["infile"] == "" ? "KPATH" : opts["infile"];

		std::ifstream ifs(infile);
		if (ifs.fail())
			throw std::runtime_error("The input file does not exist.\n");\

		interpol::kpLabelCoords kpath;
		std::string line;
		launch::Timer t;
		t.start("========= Evaluating phonon dispersion...");
		while (ifs.good())
		{
			getline(ifs, line);
			if (!line.empty())
			{
				auto line_splitted = custom_split(line, ' ');
				const std::string label = line_splitted[0].data();
				double kx = std::stod(line_splitted[1].data());
				double ky = std::stod(line_splitted[2].data());
				double kz = std::stod(line_splitted[3].data());
				
				array1D k{ kx, ky, kz };
				kpath.push_back(std::make_pair(label, k));
			}
		}
		interpol::evaluate_along_kp_path<EigenFrequencyDrawable>(kpath, bins);
		t.stop("========= Phonon dispersion is calculated. The results are written to EigenFrequency.dat.\n"
		       "          See also KLABELS for correspondance between k-point labels and k-path distance.\n");
		t.print_elapsed("\t  Phonon dispersion evaluation time: ");
		ifs.close();
	} CMD_END;

	CMD("elphmat",
		"calculates EPI matrix elecments curve along given k-path, requires input KPATH file.",
		"skies elphmat [options]:\n"
		"  Options can be given in any order. Options include:\n"
		"  Mandatory:\n"
		"  --band-ini=<num>: initial band number (starts from zero)\n"
		"  --band-fin=<num>: final band number (starts from zero)\n"
		"  --phon-branch=<num>: chosen phonon branch (starts from zero)\n"
		"  Non-mandatory:\n"
		"  --kpoint-ini=<[num, num, num]>: k-point of the initial state in crystal coordinates. Default is Gamma point [0, 0, 0]\n"
		"  --bins=<num>: number of bins between each pair of k-points in given path\n"
		"  --infile=<str>: the name of file with enumerated k-points.\n\n"
		"                  The format must be:\n"
		"                  <k-point-1 label> <num> <num> <num>\n"
		"                  <k-point-2 label> <num> <num> <num>\n"
		"                  ...................................\n"
		"                  <k-point-N label> <num> <num> <num>\n"
		"                  Here <num> after labels are crystal coordinates.\n"
		"                  Default name is 'KPATH'\n"
	) {
		if (opts["band-ini"] == "")
			throw std::runtime_error("Initial band number must be given (starts from zero)");
		if (opts["band-fin"] == "")
			throw std::runtime_error("Final band number must be given (starts from zero)");
		if (opts["phon-branch"] == "")
			throw std::runtime_error("Phonon branch number must be given (starts from zero)");

		int bins = opts["bins"] == "" ? 30 : std::stoi(opts["bins"]);
		std::string	infile = opts["infile"] == "" ? "KPATH" : opts["infile"];

		size_t nini = static_cast<size_t>(std::stoi(opts["band-ini"]));
		if (nini > EigenValue::nbands - 1)
			throw std::runtime_error("band-ini must be given in range [0, nbands)");

		size_t nfin = static_cast<size_t>(std::stoi(opts["band-fin"]));
		if (nfin > EigenValue::nbands - 1)
			throw std::runtime_error("band-fin must be given in range [0, nbands)");

		size_t nu = static_cast<size_t>(std::stoi(opts["phon-branch"]));
		if (nu > EigenFrequency::nmodes - 1)
			throw std::runtime_error("phon-branch must be given in range [0, nmodes)");

		array1D kInit{ 0.0, 0.0, 0.0 };

		std::ifstream ifs(infile);
		if (ifs.fail())
			throw std::runtime_error("The input file does not exist.\n");\

		interpol::kpLabelCoords kpath;
		std::string line;
		launch::Timer t;
		t.start("========= Evaluating EPI matrix elements along k-path...");
		while (ifs.good())
		{
			getline(ifs, line);
			if (!line.empty())
			{
				auto line_splitted = custom_split(line, ' ');
				const std::string label = line_splitted[0].data();
				double kx = std::stod(line_splitted[1].data());
				double ky = std::stod(line_splitted[2].data());
				double kz = std::stod(line_splitted[3].data());
				
				array1D k{ kx, ky, kz };
				kpath.push_back(std::make_pair(label, k));
			}
		}
		interpol::evaluate_along_kp_path<EPHMatrixDrawable>(kpath, bins, kInit, nu, nini, nfin);
		t.stop("========= EPI matrix elements along k-path are calculated. The results are written to EPHMatrix.dat.\n"
		       "          See also KLABELS for correspondance between k-point labels and k-path distance.\n");
		t.print_elapsed("\t  EPI matrix elements evaluation time: ");
		ifs.close();
	} CMD_END;

	CMD("velocs",
		"calculates velocs \"bands\" structure along given k-path",
		"skies velocs [options]:\n"
		"  Options can be given in any order. Options include:\n"
		"  --cart=<x, y or z>: cartesian component to be evaluated\n"
		"  --bins=<num>: number of bins between each pair of k-points in given path\n"
		"  --infile=<str>: the name of file with enumerated k-points.\n\n"
		"                  The format must be:\n"
		"                  <k-point-1 label> <num> <num> <num>\n"
		"                  <k-point-2 label> <num> <num> <num>\n"
		"                  ...................................\n"
		"                  <k-point-N label> <num> <num> <num>\n"
		"                  Here <num> after labels are crystal coordinates.\n"
		"                  Default name is 'KPATH'\n"
	) {
		int bins = opts["bins"] == "" ? 30 : std::stoi(opts["bins"]);
		std::string	infile = opts["infile"] == "" ? "KPATH" : opts["infile"];

		std::ifstream ifs(infile);
		if (ifs.fail())
			throw std::runtime_error("The input file does not exist.\n");\

		interpol::kpLabelCoords kpath;
		std::string line;
		launch::Timer t;
		t.start("========= Evaluating velocities \"band\" structure...");
		while (ifs.good())
		{
			getline(ifs, line);
			if (!line.empty())
			{
				auto line_splitted = custom_split(line, ' ');
				const std::string label = line_splitted[0].data();
				double kx = std::stod(line_splitted[1].data());
				double ky = std::stod(line_splitted[2].data());
				double kz = std::stod(line_splitted[3].data());
				
				array1D k{ kx, ky, kz };
				kpath.push_back(std::make_pair(label, k));
			}
		}

		char cart = opts["cart"] == "" ? 'x' : *opts["cart"].c_str();

		interpol::evaluate_along_kp_path<VelocitiesDrawable>(kpath, bins, cart);
		t.stop("========= Velocities \"band\" structure is calculated. The results are written to Velocities.dat.\n"
		       "          See also KLABELS for correspondance between k-point labels and k-path distance.\n");
		t.print_elapsed("\t  Velocities \"band\" structure evaluation time: ");
		ifs.close();
	} CMD_END;

	CMD("a2f",
		"calculates transport spectral function a2F",
		"  Options can be given in any order. Options include:\n"
		"  Mandatory:\n"
		"  --eF=<num>: Fermi level\n"
		"  --signs=[s,s']: signs in the a2f formula brackets which contain Fermi surface harmonics\n"
        "  --kgrid=<[nk1,nk2,nk3]>: Monchorst-Pack grid of k-points (nk1 x nk2 x nk3) in the 1st BZ\n"
        "  --qgrid=<[nk1,nk2,nk3]>: Monchorst-Pack grid of k-points (nk1 x nk2 x nk3) in the 1st BZ\n"
		"  Non-mandatory:\n"
        "  --epsilons=<[low,high]>: electron energy range to bracket, from low to high (in eV relative to Fermi level). Default is [-0.5,0.5]\n"
        "  --omegas=<[low,high]>: phonon frequency range to bracket, from low to high (in meV). Default is [0.1, 50.0]\n"
		"  --eps_bins=<num>: number of uniformly selected electron energies from given epsilons. Used for integration in elastic formulas. Default value is 30\n"
		"  --oms_bins=<num>: number of uniformly selected phonon energies from given omegas. Default value is 200\n"
		"  --alpha=<x, y or z>: cartesian index of electronic velocities to be used in spectral function calculation. Default is x\n"
		"  --beta=<x, y or z>: cartesian index of electronic velocities to be used in spectral function calculation. Default is x\n"
		"  --sampling=<str>: type of smearing for delta-functions approximation. Supported values are: 'gs' (gaussians), 'fd' (fermi-dirac derivative). Default value is 'fd'\n"
		"  --el_smearing=<sigma>: value of electron energies smearing in eV for chosen sampling. Default value is 0.03 eV\n"
		"  --ph_smearing=<sigma>: value of phonon frequencies smearing in eV for chosen sampling. Default value is 0.0005\n"
		"  --bands=<[low,high]>: range of band numbers (starting from 0) to consider in a2F calculation\n. Default is [0, nbnd - 1], where nbnd is total number of bands\n"
		"  --Te=<num>: electronic temperature in eV to smear transport DOS. Default value is 0.258\n"
        //"  --filename=<str>: the name of the output file. Default filename is 'SpecFunc_pp_xx.dat'\n"
        "  --continue: if specified, continue calculation of spectral function. The filename which contains an inerrupted calculation is must have name 'LambdaTr_ss'_{alpha}_{beta}.dat'\n"
		"  --tetra: if given, use doubly constrained tetrahedron method for BZ sampling\n"
	) {
		if (opts["signs"] == "")
			throw std::runtime_error("a2f command called without signs specification\n");
		auto signs = parse_vector_of_numbers<int>(opts["signs"], "signs", 2);
		int sign    = signs[0];
		int sign_pr = signs[1];

		if ((std::abs(sign) != 1) || (std::abs(sign_pr) != 1))
			throw std::runtime_error("signs must take values of +1 or -1\n");

		char alpha = opts["alpha"] == "" ? 'x' : *opts["alpha"].c_str();
		char beta  = opts["beta"]  == "" ? 'x' : *opts["beta"].c_str();

		std::string filename = "SpecFunc_";
		if (sign > 0)    filename += std::string{'p'};
		else		     filename += std::string{'m'};
		if (sign_pr > 0) filename += std::string{'p'};
		else		     filename += std::string{'m'};
		filename += std::string{'_'};
		filename += std::string{alpha};
		filename += std::string{beta};
		filename += std::string{".dat"};

		double low_freq{ 0.1 / 1000.0 };
        double high_freq{ 50.0 / 1000.0 };
		if (opts["omegas"] != "")
        {
			auto omegas = parse_vector_of_numbers<double>(opts["omegas"], "omegas", 2);
			low_freq  = omegas[0] / 1000.0;
			high_freq = omegas[1] / 1000.0;
        }

		int noms{ 200 };
		if (opts["oms_bins"] != "")
			noms = stoi(opts["oms_bins"]);
		auto omegas = create_range(low_freq, high_freq, noms);

		bool is_continue = opts["continue"] == "true" ? true : false;
		if (is_continue)
		{
			std::string cont_filename = "LambdaTr_";
			if (sign > 0)    cont_filename += std::string{'p'};
			else		     cont_filename += std::string{'m'};
			if (sign_pr > 0) cont_filename += std::string{'p'};
			else		     cont_filename += std::string{'m'};
			cont_filename += std::string{'_'};
			cont_filename += std::string{alpha};
			cont_filename += std::string{beta};
			cont_filename += std::string{".dat"};

			auto a2f = SpecFunc(cont_filename);

			launch::Timer t;
    		t.start("========= Continuation of interrupted spectral function a2F calculation... See file " + cont_filename + " for details\n");
			calc_spec_func(a2f, omegas, filename);

			t.stop("========= Transport spectral function a2F is evaluated. The results are written to " + filename);
			t.print_elapsed("\ta2F evaluation time: ");

			return cmd_found;
		}

		if (opts["eF"] == "")
			throw std::runtime_error("a2f command called without Fermi level specification\n");
		EigenValue::eF = std::stod(opts["eF"]);

		size_t low_band{ 0 };
		size_t high_band{ EigenValue::nbands - 1 };
		if (opts["bands"] != "")
		{
			auto bands = parse_vector_of_numbers<size_t>(opts["bands"], "bands", 2);
			low_band  = bands[0];
			high_band = bands[1];
		}

		if (opts["kgrid"] == "")
			throw std::runtime_error("k-point grid is not specified");
		auto kpgrid = parse_vector_of_numbers<size_t>(opts["kgrid"], "kgrid", 3);

		if (opts["qgrid"] == "")
			throw std::runtime_error("q-point grid is not specified");
		auto qpgrid = parse_vector_of_numbers<size_t>(opts["qgrid"], "qgrid", 3);

		std::string sampling_str = "fd"; 
        if (opts["sampling"] != "")
		    sampling_str = opts["sampling"];
		auto sampling = bzsampling::switch_sampling(sampling_str);
		
        double el_smearing = 0.03;
        if (opts["el_smearing"] != "")
		    el_smearing = stod(opts["el_smearing"]);

		double ph_smearing = 0.0005;
        if (opts["ph_smearing"] != "")
		    ph_smearing = stod(opts["ph_smearing"]);

		double low_energy{ -0.5 };
        double high_energy{ 0.5 };
		if (opts["epsilons"] != "")
        {
			auto epsilons = parse_vector_of_numbers<double>(opts["epsilons"], "epsilons", 2);
			low_energy  = epsilons[0];
			high_energy = epsilons[1];
        }

		int neps{ 30 };
		if (opts["eps_bins"] != "")
			neps = stoi(opts["eps_bins"]);
		auto epsilons = create_range(low_energy, high_energy, neps);

		double Te = opts["Te"] == "" ? 0.258 : stod(opts["Te"]);

		bool is_tetra = opts["tetra"] == "true" ? true : false;

		launch::Timer t;
		if (opts["epsilons"] == "")
		{
			auto a2f = SpecFunc(kpgrid,
								qpgrid,
								sampling,
								sampling,
								el_smearing,
								ph_smearing,
								sign,
								sign_pr,
								Te,
								alpha,
								beta,
								is_tetra
			);
			a2f.set_low_band(low_band);
			a2f.set_high_band(high_band);

			a2f.set_type_of_el_smear(sampling_str);
			a2f.set_type_of_ph_smear(sampling_str);

    		t.start("========= Evaluating transport spectral function a2F...");
			calc_spec_func(a2f, omegas, filename);
		}
		else
		{
			auto a2f = SpecFunc(kpgrid,
                                qpgrid,
                                epsilons,
                                sampling,
								sampling,
								el_smearing,
								ph_smearing,
                                sign,
                                sign_pr,
								Te,
								alpha,
								beta,
								is_tetra
            );
			a2f.set_low_band(low_band);
			a2f.set_high_band(high_band);

			a2f.set_type_of_el_smear(sampling_str);
			a2f.set_type_of_ph_smear(sampling_str);

			t.start("========= Evaluating transport spectral function a2F...");
			calc_spec_func(a2f, omegas, filename);
		}

		t.stop("========= Transport spectral function a2F is evaluated. The results are written to " + filename);
		t.print_elapsed("\ta2F evaluation time: ");
	} CMD_END;

	CMD("resist",
		"calculates electrical resistivity",
		"skies resist [options]:\n"
		"  Options can be given in any order. Options include:\n"
		"  Mandatory:\n"
        "  --range=<[low, high]>: electronic temperature range to calculate electrical resistivity\n"
		"  Non-mandatory:\n"
		"  --bins=<num>: how many temperature values in the given range. Default value is 200\n"
        //"  --infile=<str>: the name of the input file which contains transport calculated spectral function. Default filename is 'SpecFunc_plus.dat'\n"
        //"  --outfile=<str>: the name of the output file. Default filename is 'resist_{alpha}_{beta}.dat, where alpha and beta are cartesian components of velocities'\n"
		"  --elastic: if given calculates electrical resistivity using elastic formula in Allen's approach. Infile must contain spectral functions at least for 2 electron energies\n"
		"  --iont-range=<[low, high]>: ionic temperature range to calculate electrical resistivity\n"
		"  --iont-bins=<num>: how many temperature values in the given ion temperatures range. Default value is 200\n"
		"  --alpha=<x, y or z>: first cartesian index in resistivity tensor component. Default is 'x'\n"
		"  --beta=<x, y or z>: first cartesian index in resistivity tensor component. Default is 'x'\n"
	) {
		if (opts["range"] == "")
			throw std::runtime_error("The range of temperatures is not specified");
		array1D range = parse_vector_of_numbers<double>(opts["range"], "range", 2);
		
		auto low_temp  = range[0];
		auto high_temp = range[1];

		size_t bins = 200;
        if (opts["bins"] != "")
		    bins = static_cast<size_t>(stoi(opts["bins"]));

		auto Temps = create_range(low_temp, high_temp, bins);

		bool is_elastic = opts["elastic"] == "true" ? true : false;

		char alpha = opts["alpha"] == "" ? 'x' : *opts["alpha"].c_str();
		char beta  = opts["beta"]  == "" ? 'x' : *opts["beta"].c_str();
		std::string	a2f_fnm = "SpecFunc_pp_";
		a2f_fnm += std::string{alpha};
		a2f_fnm += std::string{beta};
		a2f_fnm += std::string{".dat"};

		std::string	outfile = "Resist_";
		outfile += std::string{alpha};
		outfile += std::string{beta};
		outfile += std::string{".dat"};

		array1D ion_Temps;
		if (opts["iont-range"] != "")
		{
			auto iont_range = parse_vector_of_numbers<double>(opts["iont-range"], "iont-range", 2);
			auto ion_low_temp  = iont_range[0];
			auto ion_high_temp = iont_range[1];
			size_t iont_bins = 200;
        	if (opts["iont-bins"] != "")
			{
		    	iont_bins = static_cast<size_t>(stoi(opts["iont-bins"]));
			}
			ion_Temps = create_range(ion_low_temp, ion_high_temp, iont_bins);
		}

		launch::Timer t;
		t.start("========= Evaluating electrical resistivity...");
		if (is_elastic)
		{
			transport::calc_elec_cond_elastic(Temps, ion_Temps, a2f_fnm.c_str(), outfile.c_str(), skies::Lattprotocol::latt_volume);
		}
		else
		{
			transport::calc_elec_cond_inelastic(Temps, a2f_fnm.c_str(), outfile.c_str(), skies::Lattprotocol::latt_volume);
		}
		t.stop("========= Electrical resistivity is evaluated. The results are written to " + outfile);
		t.print_elapsed("\t  Electrical resistivity evaluation time: ");
	} CMD_END;

	CMD("thermal-resist",
		"calculates thermal resistivity",
		"skies thermal-resist [options]:\n"
		"  Options can be given in any order. Options include:\n"
		"  Mandatory:\n"
        "  --range=<[low, high]>: temperature range to calculate thermal resistivity (in K)\n"
		"  Non-mandatory:\n"
		"  --bins=<num>: how many temperature values in the given range. Default value is 200\n"
        //"  --infile-plus=<str>: the name of the input file which contains transport calculated spectral function for +1 sign. Default filename is 'SpecFunc_plus.dat'\n"
        //"  --infile-minus=<str>: the name of the input file which contains transport calculated spectral function for -1 sign. Default filename is 'SpecFunc_minus.dat'\n"
        //"  --outfile=<str>: the name of the output file. Default filename is 'thermal-'\n"
		"  --elastic: if given calculates thermal resistivity using elastic formula in Allen's approach. Infile must contain spectral functions at least for 2 electron energies\n"
		"  --iont-range=<[low, high]>: ionic temperature range to calculate electrical resistivity\n"
		"  --iont-bins=<num>: how many temperature values in the given ion temperatures range. Default value is 200\n"
		"  --alpha=<x, y or z>: first cartesian index in resistivity tensor component. Default is 'x'\n"
		"  --beta=<x, y or z>: first cartesian index in resistivity tensor component. Default is 'x'\n"
	) {
		if (opts["range"] == "")
			throw std::runtime_error("The range of temperatures is not specified");
		array1D range = parse_vector_of_numbers<double>(opts["range"], "range", 2);
		
		auto low_temp  = range[0];
		auto high_temp = range[1];

		size_t bins = 200;
        if (opts["bins"] != "")
		    bins = static_cast<size_t>(stoi(opts["bins"]));

		auto Temps = create_range(low_temp, high_temp, bins);

		bool is_elastic = opts["elastic"] == "true" ? true : false;

		char alpha = opts["alpha"] == "" ? 'x' : *opts["alpha"].c_str();
		char beta  = opts["beta"]  == "" ? 'x' : *opts["beta"].c_str();

		std::string	a2f_plus_fnm = "SpecFunc_pp_";
		a2f_plus_fnm += std::string{alpha};
		a2f_plus_fnm += std::string{beta};
		a2f_plus_fnm += std::string{".dat"};

		std::string	a2f_minus_fnm = "SpecFunc_mm_";
		a2f_minus_fnm += std::string{alpha};
		a2f_minus_fnm += std::string{beta};
		a2f_minus_fnm += std::string{".dat"};

		std::string	outfile = "ThermalResist_";
		outfile += std::string{alpha};
		outfile += std::string{beta};
		outfile += std::string{".dat"};

		array1D ion_Temps;
		if (opts["iont-range"] != "")
		{
			auto iont_range = parse_vector_of_numbers<double>(opts["iont-range"], "iont-range", 2);
			auto ion_low_temp  = iont_range[0];
			auto ion_high_temp = iont_range[1];
			size_t iont_bins = 200;
        	if (opts["iont-bins"] != "")
			{
		    	iont_bins = static_cast<size_t>(stoi(opts["iont-bins"]));
			}
			ion_Temps = create_range(ion_low_temp, ion_high_temp, iont_bins);
		}

		bool add_odd = opts["add_odd"] == "true" ? true : false;
		std::string	a2f_pm_fnm;
		std::string	a2f_mp_fnm;
		assert(a2f_pm_fnm.size() == 0);
		assert(a2f_mp_fnm.size() == 0);
		if (add_odd)
		{
			a2f_pm_fnm = "SpecFunc_pm_";
			a2f_pm_fnm += std::string{alpha};
			a2f_pm_fnm += std::string{beta};
			a2f_pm_fnm += std::string{".dat"};

			a2f_mp_fnm = "SpecFunc_mp_";
			a2f_mp_fnm += std::string{alpha};
			a2f_mp_fnm += std::string{beta};
			a2f_mp_fnm += std::string{".dat"};
		}

		launch::Timer t;
		t.start("========= Evaluating thermal conductivity...");
		if (is_elastic)
		{
			transport::calc_therm_cond_elastic(Temps, ion_Temps, a2f_plus_fnm.c_str(), a2f_minus_fnm.c_str(),
				outfile.c_str(), skies::Lattprotocol::latt_volume, a2f_pm_fnm, a2f_mp_fnm);
		}
		else
		{
			transport::calc_therm_cond_inelastic(Temps, a2f_plus_fnm.c_str(), a2f_minus_fnm.c_str(),
				outfile.c_str(), skies::Lattprotocol::latt_volume);
		}
		t.stop("========= Thermal conductivity is evaluated. The results are written to " + outfile);
		t.print_elapsed("\tThermal conductivity evaluation time: ");
	} CMD_END;

	CMD("seebeck",
		"calculates Seebeck coefficient",
		"skies seebeck [options]:\n"
		"  Options can be given in any order. Options include:\n"
		"  Mandatory:\n"
        "  --range=<[low,high]>: temperature range to calculate thermal resistivity (in K)\n"
		"  Non-mandatory:\n"
		"  --bins=<num>: how many temperature values in the given range. Default value is 200\n"
        //"  --infile-plus=<str>: the name of the input file which contains transport calculated spectral function for +1 sign. Default filename is 'SpecFunc_plus.dat'\n"
        //"  --infile-minus=<str>: the name of the input file which contains transport calculated spectral function for -1 sign. Default filename is 'SpecFunc_minus.dat'\n"
        //"  --outfile=<str>: the name of the output file. Default filename is 'seebeck.dat'\n"
		"  --iont-range=<[low, high]>: ionic temperature range to calculate electrical resistivity\n"
		"  --iont-bins=<num>: how many temperature values in the given ion temperatures range. Default value is 200\n"
		"  --alpha=<x, y or z>: first cartesian index in resistivity tensor component. Default is 'x'\n"
		"  --beta=<x, y or z>: first cartesian index in resistivity tensor component. Default is 'x'\n"
		"  --add_odd: if given add odd transport specral function corrections\n"
	) {
		if (opts["range"] == "")
			throw std::runtime_error("The range of temperatures is not specified");
		array1D range = parse_vector_of_numbers<double>(opts["range"], "range", 2);
		
		auto low_temp  = range[0];
		auto high_temp = range[1];

		size_t bins = 200;
        if (opts["bins"] != "")
		    bins = static_cast<size_t>(stoi(opts["bins"]));

		auto Temps = create_range(low_temp, high_temp, bins);

		char alpha = opts["alpha"] == "" ? 'x' : *opts["alpha"].c_str();
		char beta  = opts["beta"]  == "" ? 'x' : *opts["beta"].c_str();

		std::string	a2f_plus_fnm = "SpecFunc_pp_";
		a2f_plus_fnm += std::string{alpha};
		a2f_plus_fnm += std::string{beta};
		a2f_plus_fnm += std::string{".dat"};

		std::string	a2f_minus_fnm = "SpecFunc_mm_";
		a2f_minus_fnm += std::string{alpha};
		a2f_minus_fnm += std::string{beta};
		a2f_minus_fnm += std::string{".dat"};

		std::string	outfile = "Seebeck_";
		outfile += std::string{alpha};
		outfile += std::string{beta};
		outfile += std::string{".dat"};

		array1D ion_Temps;
		if (opts["iont-range"] != "")
		{
			auto iont_range = parse_vector_of_numbers<double>(opts["iont-range"], "iont-range", 2);
			auto ion_low_temp  = iont_range[0];
			auto ion_high_temp = iont_range[1];
			size_t iont_bins = 200;
        	if (opts["iont-bins"] != "")
			{
		    	iont_bins = static_cast<size_t>(stoi(opts["iont-bins"]));
			}
			ion_Temps = create_range(ion_low_temp, ion_high_temp, iont_bins);
		}

		bool add_odd = opts["add_odd"] == "true" ? true : false;
		std::string	a2f_pm_fnm;
		std::string	a2f_mp_fnm;
		assert(a2f_pm_fnm.size() == 0);
		assert(a2f_mp_fnm.size() == 0);
		if (add_odd)
		{
			a2f_pm_fnm = "SpecFunc_pm_";
			a2f_pm_fnm += std::string{alpha};
			a2f_pm_fnm += std::string{beta};
			a2f_pm_fnm += std::string{".dat"};

			a2f_mp_fnm = "SpecFunc_mp_";
			a2f_mp_fnm += std::string{alpha};
			a2f_mp_fnm += std::string{beta};
			a2f_mp_fnm += std::string{".dat"};
		}
		launch::Timer t;
		t.start("========= Evaluating Seebeck coefficient...");

		transport::calc_seebeck_elastic(Temps, ion_Temps, a2f_plus_fnm.c_str(), a2f_minus_fnm.c_str(),
			outfile.c_str(), skies::Lattprotocol::latt_volume, a2f_pm_fnm, a2f_mp_fnm);

		t.stop("========= Seebeck coefficient is evaluated. The results are written to " + outfile);
		t.print_elapsed("\t  Seebeck coefficient evaluation time: ");
	} CMD_END;

    return cmd_found;
}

} // launch
} // skies