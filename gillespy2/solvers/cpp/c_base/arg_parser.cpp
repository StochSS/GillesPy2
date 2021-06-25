#include <string>
#include <sstream>
#include <iostream>

#include "template.h"
#include "arg_parser.h"

char ArgParser::match_arg(std::string &token)
{
	if (!token.compare("--timesteps"))
	{
		return 't';
	}

	if (!token.compare("--end")) 
	{
		return 'e';
	}

	if (!token.compare("--seed"))
	{
		return 's';
	}

	if (!token.compare("--switch_tol"))
	{
		return 'S';
	}

	if (!token.compare("--increment"))
	{
		return 'i';
	}

	if (!token.compare("--init_pop"))
	{
		return 'I';
	}

	if (!token.compare("--parameters"))
	{
		return 'p';
	}

	if (!token.compare("--trajectories")) 
	{
		return 'T';
	}

	if (!token.compare("--tau_tol"))
	{
		return 'l';
	}

	else
	{
		return 0;
	}
}

ArgParser::ArgParser(int argc, char *argv[])
{
	std::stringstream arg_stream;

	for (int i = 1; i < argc; i++) 
	{
		char arg_key;
		std::string token(argv[i]);

		// If the token does not start with a dash then it is not a valid argument token. Seek forward.
		if (token.rfind("-", 0) == std::string::npos) 
		{ 
			continue;
		}

		// If the arg-key starts with two dashes: '--' then it's in long form.
		if (token.rfind("--", 0) == 0)
		{
			arg_key = match_arg(token);
		}

		// Token is of short form, we may want to do some validation here later.
		else
		{
			arg_key = token[1];
		}

		// Create a stringstream for datatype parsing.
		arg_stream << argv[i + 1];
		std::stringstream arg_stream(argv[i + 1]);

		switch (arg_key) 
		{
			case 't':
				arg_stream >> timesteps;
				break;

			case 'e':
				arg_stream >> end;
				break;

			case 's':
				arg_stream >> seed;
				break;

			case 'S':
				arg_stream >> switch_tol;
				break;

			case 'i':
				arg_stream >> increment;
				break;

			case 'I':
				Gillespy::map_variable_populations(arg_stream);
				break;

			case 'p':
				Gillespy::map_variable_parameters(arg_stream);
				break;

			case 'T':
				arg_stream >> trajectories;
				break;

			case 'l':
				arg_stream >> tau_tol;
				break;

			default:
				std::cerr << usage << std::endl;
				exit(0);
		}
	}

	// Negative or 0 increment is invalid.
	// If such an increment is provided, compute a valid increment manually.
	increment = increment > 0 ? increment : end / (timesteps - 1);
}

ArgParser::~ArgParser()
{
}