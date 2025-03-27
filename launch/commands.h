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

#include <launch/types.h>

namespace skies { namespace launch {

/**
 * @brief Enumeration representing supported command-line commands.
 *
 * Each value corresponds to a specific command that the program can execute.
 */
enum class CMD {
	dummy,          ///< Placeholder for an invalid or unspecified command.
	help,           ///< Displays general help information.
	list,           ///< Lists all available commands.
	dos,            ///< Computes density of states (DOS).
	phdos,          ///< Computes phonon density of states (phDOS).
	trdos,          ///< Computes transport density of states (trDOS).
	bands,          ///< Computes electronic band structures.
	phonons,        ///< Computes phonon properties.
	elphmat,        ///< Computes electron-phonon coupling matrices.
	velocs,         ///< Computes electronic velocities.
	a2f,            ///< Computes electron-phonon spectral function (a2F).
	resist,         ///< Computes electrical resistivity.
	thermal_cond    ///< Computes thermal conductivity.
};

/**
 * @brief Converts a string representation of a command to its corresponding `CMD` enumeration value.
 *
 * @param cmd String representation of the command (e.g., "dos", "phdos").
 * @return CMD The corresponding `CMD` enumeration value.
 */
CMD str_to_CMD(const std::string& cmd);

/**
 * @brief Converts a `CMD` enumeration value to its string representation.
 *
 * @param cmd The `CMD` enumeration value to convert.
 * @return std::string The string representation of the command.
 */
std::string CMD_to_str(CMD cmd);

/**
 * @brief Displays help information for a specific command.
 *
 * @param cmd The command for which help information is displayed.
 */
void help_for_cmd(CMD cmd);

/**
 * @brief Resolves and executes a command based on the provided arguments and options.
 *
 * @param cmd The command to execute.
 * @param args Parsed command-line arguments.
 * @param opts Parsed command-line options.
 */
void resolve_cmd(CMD cmd, const TArgs& args, TOpts& opts);

/**
 * @brief Displays help information for the `dos` command.
 */
void help_for_dos();

/**
 * @brief Displays help information for the `phdos` command.
 */
void help_for_phdos();

/**
 * @brief Displays help information for the `trdos` command.
 */
void help_for_trdos();

/**
 * @brief Displays help information for the `bands` command.
 */
void help_for_bands();

/**
 * @brief Displays help information for the `phonons` command.
 */
void help_for_phonons();

/**
 * @brief Displays help information for the `elphmat` command.
 */
void help_for_elphmat();

/**
 * @brief Displays help information for the `velocs` command.
 */
void help_for_velocs();

/**
 * @brief Displays help information for the `a2f` command.
 */
void help_for_a2f();

/**
 * @brief Displays help information for the `resist` command.
 */
void help_for_resist();

/**
 * @brief Displays help information for the `thermal_cond` command.
 */
void help_for_thermal_cond();

/**
 * @brief Executes the `list` command, listing all available commands.
 */
void cmd_list();

/**
 * @brief Executes the `dos` command with the provided options.
 *
 * @param opts Parsed command-line options.
 */
void cmd_dos(TOpts& opts);

/**
 * @brief Executes the `phdos` command with the provided options.
 *
 * @param opts Parsed command-line options.
 */
void cmd_phdos(TOpts& opts);

/**
 * @brief Executes the `trdos` command with the provided options.
 *
 * @param opts Parsed command-line options.
 */
void cmd_trdos(TOpts& opts);

/**
 * @brief Executes the `bands` command with the provided options.
 *
 * @param opts Parsed command-line options.
 */
void cmd_bands(TOpts& opts);

/**
 * @brief Executes the `phonons` command with the provided options.
 *
 * @param opts Parsed command-line options.
 */
void cmd_phonons(TOpts& opts);

/**
 * @brief Executes the `elphmat` command with the provided options.
 *
 * @param opts Parsed command-line options.
 */
void cmd_elphmat(TOpts& opts);

/**
 * @brief Executes the `velocs` command with the provided options.
 *
 * @param opts Parsed command-line options.
 */
void cmd_velocs(TOpts& opts);

/**
 * @brief Executes the `a2f` command with the provided options.
 *
 * @param opts Parsed command-line options.
 */
void cmd_a2f(TOpts& opts);

/**
 * @brief Executes the `resist` command with the provided options.
 *
 * @param opts Parsed command-line options.
 */
void cmd_resist(TOpts& opts);

/**
 * @brief Executes the `thermal_cond` command with the provided options.
 *
 * @param opts Parsed command-line options.
 */
void cmd_thermal_cond(TOpts& opts);

} // namespace launch
} // namespace skies