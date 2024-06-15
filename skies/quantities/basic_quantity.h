#pragma once

#include <skies/common/units.h>
#include <skies/sampling/sampling.h>
#include <skies/common/ndimarrays.h>

namespace skies { namespace quantities {

class AnyQuantity {
public:
    virtual std::string name() const = 0;
    virtual arrays::array1D interpolate_at(const arrays::array1D& q) const = 0;
    virtual ~AnyQuantity() {}
};

} // quantities
} // skies
