#pragma once

#include <launch/types.h>

namespace skies { namespace launch {

enum class CMD {
	dummy,
	help,
	list,
	dos,
	phdos,
	trdos,
	bands,
	phonons,
	elphmat,
	velocs,
	a2f,
	resist,
	thermal_cond
};

CMD str_to_CMD(const std::string& cmd);

std::string CMD_to_str(CMD cmd);

void help_for_cmd(CMD cmd);

void resolve_cmd(CMD cmd, const TArgs& args, TOpts& opts);

void help_for_dos();
void help_for_phdos();
void help_for_trdos();
void help_for_bands();
void help_for_phonons();
void help_for_elphmat();
void help_for_velocs();
void help_for_a2f();
void help_for_resist();
void help_for_thermal_cond();

void cmd_list();
void cmd_dos(TOpts& opts);
void cmd_phdos(TOpts& opts);
void cmd_trdos(TOpts& opts);
void cmd_bands(TOpts& opts);
void cmd_phonons(TOpts& opts);
void cmd_elphmat(TOpts& opts);
void cmd_velocs(TOpts& opts);
void cmd_a2f(TOpts& opts);
void cmd_resist(TOpts& opts);
void cmd_thermal_cond(TOpts& opts);

} // launch
} // skies