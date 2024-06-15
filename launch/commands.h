#pragma once

#include <vector>
#include <string>
#include <unordered_map>

namespace skies { namespace launch {

#define CMD(cmd_name, descr, how_to) \
	if(cmd == "list" && !rank) \
		std::cout << "    "  << cmd_name << ": " << descr << '\n'; \
	else if((args.size() == 1) && (cmd == "help") && (args[0] == cmd_name) && !rank) \
		{std::cout << "Usage: \n" << how_to; cmd_found = true; } \
	else if (cmd == cmd_name) { cmd_found = true;

#define CMD_END }

void exec_cmd(const std::string& cmd,
	         std::vector<std::string>& args,
	         std::unordered_map<std::string, std::string>& opts);

bool resolve_cmd(const std::string& cmd,
                 std::vector<std::string>& args,
                 std::unordered_map<std::string, std::string>& opts);

} // launch
} // skies