#pragma once

#include <vector>
#include <string>
#include <unordered_map>

namespace skies { namespace launch {

void parse_opts(int argc, char* argv[],
                std::vector<std::string>& args,
	            std::unordered_map<std::string, std::string>& opts);

} // launch
} // skies