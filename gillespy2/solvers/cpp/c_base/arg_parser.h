#include <map>
#include <string>
#include <vector>
#include <iostream>
// #include <getopt.h>
#include <string.h>
// #include <stdlib.h>

#include "template.h"

class ArgParser
{
private:
    char* usage = "\
usage: [simulation.out] \
[-t|--timesteps] <int>\
[-e|--end] <int|double> \
[-s|--seed] <int> \
[-S|--switch_tol] <double> \
[-i|--increment] <int|double> \
[-I|--init_pop] <int>... \
[-p|--parameters] <int|double>... \
[-T|--trajectories] <int>... \
[--tau_tol] <double> \n\
";
    char match_arg(char* opt)
    {
        if (!strcmp(opt, "--timesteps"))
            return 't';
        if (!strcmp(opt, "--end"))
            return 'e';
        if (!strcmp(opt, "--seed"))
            return 's';
        if (!strcmp(opt, "--switch_tol"))
            return 'S';
        if (!strcmp(opt, "--increment")) //double
            return 'i';
        if (!strcmp(opt, "--init_pop"))
            return 'I';
        if (!strcmp(opt, "--parameters"))
            return 'p';
        if (!strcmp(opt, "--trajectories"))
            return 'T';
        if (!strcmp(opt, "--tau_tol"))
            return 'l';
        else{
            printf(usage);
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

    ArgParser(const int argc, char *argv[]);
    ~ArgParser();
};

ArgParser::ArgParser(int argc, char* argv[])    
{
    int opt;

    for (int i = 1; i < argc; ++i){
        if (strlen(argv[i]) > 1 && argv[i][0] == '-' && argv[i][1] != '-') {
            opt = argv[i][1];
        }
        else if (strlen(argv[i]) > 1 && argv[i][0] == '-' && argv[i][1] == '-')
        {
            opt = match_arg(argv[i]);
        }
        else {
            continue;
        }
        std :: stringstream arg_stream(argv[i+1]);
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
                break;
            case 'l':
                arg_stream >> tau_tol;
                break;
            default:
                printf(usage);
                break;
        }        
    }

}

ArgParser::~ArgParser(){
}
