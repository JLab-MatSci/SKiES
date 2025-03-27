/*-----------------------------------------------------------------------
    * SKiES - Solver of Kinetic Equation for Solids
    * 
    * (C) 2025 Galtsov Ilya, Fokin Vladimir, Minakov Dmitry, Levashov Pavel (JIHT RAS)
    *
    * SKiES may only be utilized for non-profit research.
    * Citing appropriate sources is required when using SKiES.
    * 
    * Distribution of this file is permitted by the GNU General Public License.
    * Examine the `LICENSE' file located in the current distribution's root directory.
------------------------------------------------------------------------- */
#pragma once

#include <skies/common/units.h>
#include <skies/sampling/sampling.h>
#include <skies/common/ndimarrays.h>

namespace skies { namespace quantities {

/**
 * @class AnyQuantity
 * @brief Abstract base class for physical quantities.
 * 
 * This class provides common interfaces for all physical quantities,
 * including methods for interpolation and naming.
 */
class AnyQuantity {
public:
    virtual std::string name() const = 0;
    virtual arrays::array1D interpolate_at(const arrays::array1D& q) const = 0;
    virtual ~AnyQuantity() {}
};

} // quantities
} // skies
