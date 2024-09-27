#pragma once

#include <vector>
#include <string>
#include <unordered_map>


namespace skies { namespace launch {

namespace {
using TArgs = std::vector<std::string>;
using TOpts = std::unordered_map<std::string, std::string>;
}

} // launch
} // skies