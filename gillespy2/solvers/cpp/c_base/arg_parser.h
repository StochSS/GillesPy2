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

    double end = 0.0;
    double increment = 0.0;
    double switch_tol = 0.0;
    double tau_tol = 0.0;

    ArgParser(int argc, char *argv[]);
    ~ArgParser();
};