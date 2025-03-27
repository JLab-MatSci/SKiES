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

 #include <map>
 #include <string>
 #include <functional>
 #include <stdexcept>
 
 #include <skies/common/ndimarrays.h>
 
 namespace skies { namespace bzsampling {
 
 /**
  * @brief Alias for a sampling function type.
  *
  * A sampling function takes two `double` arguments (e.g., energy and smearing width) and returns a `double`.
  */
 using SamplingFunc = std::function<double(double, double)>;
 
 /**
  * @brief Alias for a map of sampling parameters.
  *
  * A `SamplingParams` object is a map where keys and values are strings, representing parameter names and their values.
  */
 using SamplingParams = std::map<const std::string, const std::string>;
 
 /**
  * @brief Enumeration for sampling types.
  *
  * Defines the supported types of sampling functions:
  * - `gs`: Gaussian sampling.
  * - `fd`: Fermi-Dirac sampling.
  */
 enum class SamplType {
     gs, ///< Represents Gaussian sampling.
     fd  ///< Represents Fermi-Dirac sampling.
 };
 
 /**
  * @brief Converts a string representation of a sampling type to its corresponding `SamplType` enum value.
  *
  * @param type String representation of the sampling type ("gs" or "fd").
  * @return SamplType The corresponding `SamplType` enum value.
  * @throws std::runtime_error If the input string does not match any known sampling type.
  */
 inline SamplType hash_type(const std::string& type)
 {
     if (type == "gs") return SamplType::gs;
     if (type == "fd") return SamplType::fd;
     throw std::runtime_error("Unknown type of sampling. Please enter one of 'gs' or 'fd'");
 }
 
 /**
  * @brief Computes the Fermi-Dirac distribution function.
  *
  * @param x Energy value.
  * @param sigma Smearing width.
  * @return double The Fermi-Dirac distribution value.
  */
 double fermi_dirac(double x, double sigma);
 
 /**
  * @brief Computes the Bose-Einstein distribution function.
  *
  * @param x Energy value.
  * @param sigma Smearing width.
  * @return double The Bose-Einstein distribution value.
  */
 double bose_einstein(double x, double sigma);
 
 /**
  * @brief Computes the derivative of the Fermi-Dirac distribution function.
  *
  * @param x Energy value.
  * @param sigma Smearing width.
  * @return double The derivative of the Fermi-Dirac distribution value.
  */
 double fd_derivative(double x, double sigma);
 
 /**
  * @brief Computes the Gaussian distribution function.
  *
  * @param x Energy value.
  * @param sigma Smearing width.
  * @return double The Gaussian distribution value.
  */
 double gauss(double x, double sigma);
 
 /**
  * @brief Returns a sampling function based on the specified type.
  *
  * This function maps a string representation of a sampling type to its corresponding function.
  *
  * @param type String representation of the sampling type ("gs" or "fd").
  * @return SamplingFunc The corresponding sampling function.
  */
 inline SamplingFunc switch_sampling(std::string type)
 {
     switch (hash_type(type))
     {
         case SamplType::gs:
             return gauss;
         case SamplType::fd:
             return fd_derivative;
         default:
             break;
     }
     return SamplingFunc{};
 }
 
 /**
  * @brief Smears a quantity using the Fermi-Dirac distribution.
  *
  * This function applies Fermi-Dirac smearing to a given quantity over a range of values.
  *
  * @param quan Array of quantities to be smeared.
  * @param range Array of range values (e.g., energy grid).
  * @param e Reference energy value.
  * @param sigma Smearing width.
  * @return double The smeared value.
  */
 double smear_with_fd(const arrays::array1D& quan, const arrays::array1D& range, double e, double sigma);
 
 } // namespace bzsampling
 } // namespace skies