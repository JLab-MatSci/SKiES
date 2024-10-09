#pragma once

#include <launch/types.h>

namespace skies { namespace launch {

void parse_opts(int argc, char* argv[], TArgs& args, TOpts& opts);
void check_opt_length(const std::string& opt, bool& wait_for_opts, TOpts& opts);

} // launch
} // skies