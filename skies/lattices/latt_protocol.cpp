#include <skies/lattices/latt_protocol.h>
#include <skies/common/units.h>

using namespace skies::arrays;

using LP = skies::Lattprotocol;

double      LP::latt_const;
array2D     LP::latt_coords;
double      LP::latt_volume;
LP::Lattice LP::latt_type;

#ifdef __cplusplus
extern "C" {
#endif

void fillLattice(double* lat_const, double* unit_cell_vol, double coords[3][3]) {
    LP::latt_const = *lat_const * skies::units::bohr_in_A; // in [A]
    LP::latt_coords.resize(3, array1D(3));
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            LP::latt_coords[j][i] = LP::latt_const * coords[i][j];
    LP::latt_volume = *unit_cell_vol; // in [bohr^3]
}

#ifdef __cplusplus
}
#endif

namespace skies {
Lattprotocol::Lattprotocol(Lattice lat, double a)
{
    latt_const = a;
    latt_type = lat;
    switch (lat)
    {
    case Lattice::fcc:
        latt_coords = {
            { a / 2., 0.0000, a / 2. },
            { a / 2., a / 2., 0.0000 },
            { 0.0000, a / 2., a / 2. }
        };
        break;

    case Lattice::bcc:
        latt_coords = {
            { a / 2., - a / 2., a / 2.},
            { a / 2.,   a / 2., a / 2.},
            {- a / 2.,  a / 2., a / 2.}
        };
        break;

    case Lattice::sc:
        latt_coords = {
            { a, 0.0000, 0.0000 },
            { 0.0000, a, 0.0000 },
            { 0.0000, 0.0000, a }
        };
        break; 

    default:
        break;
    }
}

array2D Lattprotocol::calc_inv_cell()
{
    auto inv = calc_inv_3d(Lattprotocol::latt_coords);
    for (auto&& v : inv)
        for (auto&& x : v)
            x *= 2.0 * atan(1.0) * 4.0;
    return transpose(inv);
}

} // skies
