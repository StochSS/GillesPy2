/*
 * GillesPy2 is a modeling toolkit for biochemical simulation.
 * Copyright (C) 2019-2021 GillesPy2 developers.
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

#pragma once

class ArgParser
{
private:
    std::string usage = "\
        usage: [simulation.out] \n\
        [-t|--timesteps] <int> \n\
        [-e|--end] <int|double> \n\
        [-s|--seed] <int> \n\
        [-S|--switch_tol] <double> \n\
        [-i|--increment] <int|double> \n\
        [-I|--init_pop] <int>... \n\
        [-p|--parameters] <int|double>... \n\
        [-T|--trajectories] <int>... \n\
        [-l|--tau_tol] <double> \n\
        ";

    char match_arg(std::string &token);

public:
    int seed = -1;
    int timesteps = 0;
    int trajectories = 0;
    int output_interval = 1;

    double end = 0.0;
    double increment = 0.0;
    double switch_tol = 0.0;
    double tau_tol = 0.03;

    ArgParser(int argc, char *argv[]);
    ~ArgParser();
};