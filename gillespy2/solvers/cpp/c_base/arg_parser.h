#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <string.h>

#include "template.h"

class ArgParser
{
private:
    std::string usage = "\
        usage: [simulation.out] \
        [-t|--timesteps] <int>\
        [-e|--end] <int|double> \
        [-s|--seed] <int> \
        [-S|--switch_tol] <double> \
        [-i|--increment] <int|double> \
        [-I|--init_pop] <int>... \
        [-p|--parameters] <int|double>... \
        [-T|--trajectories] <int>... \
        [-l|--tau_tol] <double> \n\
        ";

    char match_arg(std::string &opt)
    {
        if (!opt.compare("--timesteps"))
            return 't';
        if (!opt.compare("--end"))
            return 'e';
        if (!opt.compare("--seed"))
            return 's';
        if (!opt.compare("--switch_tol"))
            return 'S';
        if (!opt.compare("--increment")) //double
            return 'i';
        if (!opt.compare("--init_pop"))
            return 'I';
        if (!opt.compare("--parameters"))
            return 'p';
        if (!opt.compare("--trajectories"))
            // printf(opt.c_str());
            // fflush(0);
            return 'T';
        if (!opt.compare("--tau_tol"))
            return 'l';
        else{
            std::cout << usage << std::endl;
        }
              
    };
public:
    int trajectories = 0;
    int timesteps = 0;
    double end = 0.0;
    int seed = -1;
    double increment = 0.0;
    double switch_tol = 0.0;
    double tau_tol = 0.0;

    ArgParser(std::stringstream &arg_stream);
    ~ArgParser();
};

ArgParser::ArgParser(std::stringstream &arg_stream)    
{
    char opt;
    // printf("hello");
    // fflush(0);

    // std :: stringstream arg_stream(argv[1]);
    std::string token;
    // printf(arg_stream.str());
    // for (int i = 1; i < argc; ++i){
    while (arg_stream >> token) {
        // printf("hey");
        std::cout << token.length() << std::endl;
        // fflush(0);
        // printf(token.c_str());
        // fflush(0);
        if (token.length() > 1 && token.at(0) == '-' && token.at(1) != '-')
        {
            opt = token.at(1);
        }
        else if (token.length() > 1 && token.at(0) == '-' && token.at(1) == '-')
        {
            
            opt = match_arg(token);
        }
        else {
            continue;
        }
        switch(opt) {
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
                std::cout << trajectories << std::endl;

                fflush(0);
                break;
            case 'l':
                arg_stream >> tau_tol;
                break;
            default:
                std::cout << usage << std::endl;
                break;
        }        
    }

}

ArgParser::~ArgParser(){
}
