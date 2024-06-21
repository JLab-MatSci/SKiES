#include <skies/sampling/sampling.h>
#include <skies/lattices/latt_protocol.h>
#include <skies/quantities/elvelocs.h>

using namespace skies::arrays;

int skies::quantities::Velocities::dir_;
int skies::quantities::VelocitiesDrawable::dir_;

namespace skies { namespace quantities {

std::string VelocitiesDrawable::name() const { return "Velocities"; }

VelocitiesDrawable::VelocitiesDrawable() {}

VelocitiesDrawable::VelocitiesDrawable(char cart)
{
    auto dir_it = avail_directions.find(cart);
    if (dir_it == avail_directions.end())
        throw std::runtime_error("Wrong cartesian index given\n");
    dir_ = dir_it->second;
}

VelocitiesDrawable::~VelocitiesDrawable() {}

arrays::array1D VelocitiesDrawable::interpolate_at(const arrays::array1D& k) const
{
    assert(k.size() == 3);

    if (!EigenValue::nbands)
        throw std::runtime_error("Please call dos() or band_structre() functions first to initalize number of bands.");
    double rvels[EigenValue::nbands];
    interpVelocityAt(&k[0], &k[1], &k[2], rvels, &dir_);
    auto vels = array1D(rvels, rvels + EigenValue::nbands);
    return vels;
}

Velocities::Velocities(char cart)
{
    auto dir_it = avail_directions.find(cart);
    if (dir_it == avail_directions.end())
        throw std::runtime_error("Wrong cartesian index given\n");
    dir_ = dir_it->second;
}

Velocities::~Velocities() {}

arrays::array1D Velocities::interpolate_at(const arrays::array1D& k)
{
    assert(k.size() == 3);

    if (!EigenValue::nbands)
        throw std::runtime_error("Please call dos() or band_structre() functions first to initalize number of bands.");
    double rvels[EigenValue::nbands];
    interpVelocityAt(&k[0], &k[1], &k[2], rvels, &dir_);
    auto vels = array1D(rvels, rvels + EigenValue::nbands);
    return vels;
}

} // quantities
} // skies
