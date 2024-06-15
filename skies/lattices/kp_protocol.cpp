#include <cmath>
#include <algorithm>
#include <float.h>
#include <iostream>
#include <numeric>

#ifdef SKIES_MPI
#include <mpi.h>
#endif

#include <skies/common/alg.h>
#include <skies/lattices/kp_protocol.h>

using namespace skies::arrays;

namespace skies {

KPprotocol::KPprotocol(int n1, int n2, int n3)
    : nkpt(n1 * n2 * n3), grid(nkpt, array1D(3))
{
    array1D v1(n1);
    array1D v2(n2);
    array1D v3(n3);
    std::iota(v1.begin(), v1.end(), 1);
    std::iota(v2.begin(), v2.end(), 1);
    std::iota(v3.begin(), v3.end(), 1);

    std::for_each(v1.begin(), v1.end(), [n1] (double& i) { i = (2.0 * i - n1 - 1) / (2.0 * n1); return i; });
    std::for_each(v2.begin(), v2.end(), [n2] (double& i) { i = (2.0 * i - n2 - 1) / (2.0 * n2); return i; });
    std::for_each(v3.begin(), v3.end(), [n3] (double& i) { i = (2.0 * i - n3 - 1) / (2.0 * n3); return i; });

    int cnt{ 0 };
    for (auto i1 : v1)
        for (auto i2 : v2)
            for (auto i3 : v3)
            {
                grid[cnt][0] = i1;
                grid[cnt][1] = i2;
                grid[cnt][2] = i3;
                cnt++;
            }
}

} // skies
