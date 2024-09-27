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
        {
			if (opt.length() == 1) 
                throw std::runtime_error("Unsupported option encountered: '-'\n");
			if (opt[1] != '-')
                throw std::runtime_error("Options which start with '--' are only supported\n");
			if (opt.length() == 2)
				wait_for_opts = false;
            else
            {
				size_t p = opt.find_first_of("=");
				if (p == std::string::npos)
					opts[opt.substr(2)] = "true";
				else
					opts[opt.substr(2, p - 2)] = opt.substr(p + 1);
			}
		}
        else
		    args.push_back(argv[i]);
	}
}

} // launch
} // skies