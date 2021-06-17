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
[-l|--tau_tol] <double> \n\
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
    char opt;

    std :: stringstream arg_stream(argv[1]);
    char * token;
    // for (int i = 1; i < argc; ++i){
    while (arg_stream >> token) {
        printf("%s\n", token);
        if (strlen(token) > 1 && token[0] == '-' && token[1] != '-')
        {
            opt = token[1];
        }
        else if (strlen(token) > 1 && token[0] == '-' && token[1] == '-')
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
