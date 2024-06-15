#include <stdexcept>

#include <launch/parser.h>

namespace skies { namespace launch {

void parse_opts(int argc, char* argv[],
                  std::vector<std::string>& args,
	              std::unordered_map<std::string, std::string>& opts)
{
	bool wait_for_opts{ true };
	for (int i = 1; i < argc; ++i) 
    {
		std::string opt{ argv[i] };
		if (opt[0] == '-' && wait_for_opts)
        {
			if (opt.length() == 1) 
                throw std::runtime_error("Unsupported option encountered: '-'\n");
			if (opt[1] != '-')
                throw std::runtime_error("Options which start with '--' are only supported\n");
			if (opt.length() == 2)
            {
				// this is '--'. Switch off expect_options
				wait_for_opts = false;
			}
            else
            {
				size_t eq_pos = opt.find_first_of("=");
				if (eq_pos == std::string::npos)
					opts[opt.substr(2)] = "true";
				else
					opts[opt.substr(2, eq_pos - 2)] = opt.substr(eq_pos + 1);
			}
		}
        else
		    args.push_back(argv[i]);
	}
}

} // launch
} // skies