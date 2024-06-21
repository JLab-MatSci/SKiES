# pragma once

//#include <functional>
#include <map>

#include <skies/quantities/eigenvals.h>

#ifdef __cplusplus
extern "C" {
#endif

void interpVelocityAt(const double* kx, const double* ky, const double* kz, double[], const int* dir);

#ifdef __cplusplus
}
#endif

namespace {
    std::map<char, int> avail_directions{
        {'x', 1},
        {'y', 2},
        {'z', 3}
    };
}

namespace skies { namespace quantities {

class VelocitiesDrawable : public AnyQuantity {
public:
    virtual std::string name() const override;
    virtual arrays::array1D interpolate_at(const arrays::array1D& q) const override;

    VelocitiesDrawable();
    VelocitiesDrawable(char cart);

    virtual ~VelocitiesDrawable();

private:
    static int dir_;
};

class Velocities {
public:
    arrays::array1D interpolate_at(const arrays::array1D& q);

    Velocities(char cart);

    virtual ~Velocities();

private:
    static int dir_;
};

} // quantities
} // skies
