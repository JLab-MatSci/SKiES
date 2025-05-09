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
#include <stdexcept>

#include <launch/parser.h>

namespace skies { namespace launch {

void parse_opts(int argc, char* argv[], TArgs& args, TOpts& opts)
{
	bool wait_for_opts{ true };
	for (int i = 1; i < argc; ++i) 
    {
		std::string opt{ argv[i] };
		if (opt[0] == '-' && wait_for_opts)
			check_opt_length(opt, wait_for_opts, opts);
        else
		    args.push_back(argv[i]);
	}
}

void check_opt_length(const std::string& opt, bool& wait_for_opts, TOpts& opts)
{
	if (opt.length() == 1)
        throw std::runtime_error("Unsupported option encountered: '-'\n");
	if (opt[1] != '-')
		throw std::runtime_error("Options which start with '--' are only supported\n");
	if (opt.length() == 2)
		wait_for_opts = false;
	else
	{
		auto p = opt.find_first_of("=");
		if (p == std::string::npos) opts[opt.substr(2)] = "true";
		else opts[opt.substr(2, p - 2)] = opt.substr(p + 1);
	}
}

} // launch
} // skies