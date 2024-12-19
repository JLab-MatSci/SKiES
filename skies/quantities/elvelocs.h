/*-----------------------------------------------------------------------
    * SKiES - Solver of Kinetic Equation for Solids
    * Version 1.0.0
    * 
    * (C) 2024 Galtsov Ilya, Fokin Vladimir, Minakov Dmitry, Levashov Pavel (JIHT RAS)
    *
    * SKiES may only be utilized for non-profit research.
    * Citing appropriate sources is required when using SKiES.
    * 
    * Project homepage: https://github.com/i-Galts/SKiES
    * 
    * This file: Electron velocities calculations implementation
    * Implements interpolation of electron velocities at arbitrary q-points
    * and provides C-interface for external velocity calculations.
    * 
    * Distribution of this file is permitted by the GNU General Public License.
    * Examine the `LICENSE' file located in the current distribution's root directory.
------------------------------------------------------------------------- */
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

/**
 * @class VelocitiesDrawable
 * @brief Class for handling drawable electron velocity data.
 * 
 * This class extends the AnyQuantity class and provides methods for
 * interpolating electron velocity data at specified points.
 */
class VelocitiesDrawable : public AnyQuantity {
public:
    /**
     * @brief Returns the name of the quantity.
     * @return A string representing the name of the quantity.
     */
    std::string name() const override;

    /**
     * @brief Interpolates electron velocities at the specified q-point.
     * @param q The q-point at which to interpolate.
     * @return An array of interpolated electron velocities.
     */
    arrays::array1D interpolate_at(const arrays::array1D& q) const override;

    /**
     * @brief Default constructor for VelocitiesDrawable.
     */
    VelocitiesDrawable();

    /**
     * @brief Constructor for VelocitiesDrawable with specified cartesian direction.
     * @param cart The cartesian direction for the velocities.
     */
    VelocitiesDrawable(char cart);

    /**
     * @brief Destructor for VelocitiesDrawable.
     */
    virtual ~VelocitiesDrawable();

private:
    static int dir_; ///< Direction of the velocities
};

/**
 * @class Velocities
 * @brief Class for handling electron velocity calculations.
 * 
 * This class provides methods for interpolating electron velocities
 * at specified points.
 */
class Velocities {
public:
    /**
     * @brief Interpolates electron velocities at the specified q-point.
     * @param q The q-point at which to interpolate.
     * @return An array of interpolated electron velocities.
     */
    arrays::array1D interpolate_at(const arrays::array1D& q);

    /**
     * @brief Constructor for Velocities with specified cartesian direction.
     * @param cart The cartesian direction for the velocities.
     */
    Velocities(char cart);

    /**
     * @brief Destructor for Velocities.
     */
    virtual ~Velocities();

private:
    static int dir_; ///< Direction of the velocities
};

} // quantities
} // skies
