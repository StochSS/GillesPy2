/*
 * GillesPy2 is a modeling toolkit for biochemical simulation.
 * Copyright (C) 2019-2022 GillesPy2 developers.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <sstream>
#include <iostream>

#include "template.h"
#include "arg_parser.h"

char ArgParser::match_arg(std::string &token)
{
    if (!token.compare("--verbose"))
    {
        return 'v';
    }

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

    if (!token.compare("--interval"))
    {
        return 'V';
    }

    if (!token.compare("--rtol"))
    {
        return 'R';
    }

    if (!token.compare("--atol"))
    {
        return 'A';
    }

    if (!token.compare("--max_step"))
    {
        return 'M';
    }

    if (!token.compare("--use_root_finding"))
    {
        return 'u';
    }

    return 0;
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

        switch (arg_key) 
        {
            case 't':
                std::stringstream(argv[i + 1]) >> timesteps;
                break;

            case 'e':
                std::stringstream(argv[i + 1]) >> end;
                break;

            case 's':
                std::stringstream(argv[i + 1]) >> seed;
                break;

            case 'S':
                std::stringstream(argv[i + 1]) >> switch_tol;
                break;

            case 'i':
                std::stringstream(argv[i + 1]) >> increment;
                break;

            case 'I':
                arg_stream = std::stringstream(argv[i + 1]);
                Gillespy::map_variable_populations(arg_stream);
                break;

            case 'p':
                arg_stream = std::stringstream(argv[i + 1]);
                Gillespy::map_variable_parameters(arg_stream);
                break;

            case 'T':
                std::stringstream(argv[i + 1]) >> trajectories;
                break;

            case 'l':
                std::stringstream(argv[i + 1]) >> tau_tol;
                break;

            case 'V':
                std::stringstream(argv[i + 1]) >> output_interval;
                break;

            case 'v':
                verbose = true;
                break;

            case 'u':
                use_root_finding = true;
                break;

            case 'R':
                std::stringstream(argv[i + 1]) >> rtol;
                break;

            case 'A':
                std::stringstream(argv[i + 1]) >> atol;
                break;

            case 'M':
                std::stringstream(argv[i + 1]) >> max_step;
                break;

            default:
                std::cerr << usage << std::endl;
                exit(0);
        }
    }

    // Negative or 0 increment is invalid.
    // If such an increment is provided, compute a valid increment manually.
    increment = increment > 0 ? increment : end / (timesteps - 1);
    // Output interval must lie within the range (0, num_timesteps]
    output_interval = std::max(1, std::min(timesteps, output_interval));
}

ArgParser::~ArgParser()
{
}
