#pragma once

#ifdef SKIES_TBB
#include <execution>
#define PAR std::execution::par,
#else
#define PAR
#endif