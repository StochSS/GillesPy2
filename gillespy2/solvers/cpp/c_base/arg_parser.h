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
        if (!strcmp(opt, "--trajectories")){
            // printf("match");
            return 'T';
        }
        
        if (!strcmp(opt, "tau_tol"))
            return 0;
        else
        {
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
    // std::vector<int> init_pop = std::vector<int>();
    // std::vector<double> parameters = std::vector<double>();

    ArgParser(const int argc, char *argv[]);
    ~ArgParser();
};

ArgParser::ArgParser(int argc, char* argv[])    
{
    int opt;
    // printf("hello2\n");
    // std::cerr << "SUP";
    for (int i = 1; i < argc; ++i){
        // printf("%s strlen(argv[%d]): %d\n",argv[i],i,strlen(argv[i]));
        if (strlen(argv[i]) > 1 && argv[i][0] == '-' && argv[i][1] != '-') {
            opt = argv[i][1];
        }
        else if (strlen(argv[i]) > 1 && argv[i][0] == '-' && argv[i][1] == '-')
        {
            // printf("long\n");
            opt = match_arg(argv[i]);
        }
        else {
            continue;
        }
        // printf("%c\n", opt);
        std :: stringstream arg_stream(argv[i+1]);
        switch(opt) {
            case 't':
                arg_stream >> timesteps;
                // printf("timesteps: %d\n", timesteps);
                break;
            case 'e':
                arg_stream >> end;
                // printf("end: %f\n", end);
                // end = strtod(argv[i+1],0);
                break;
            case 's':
                arg_stream >> seed;
                // seed = atoi(argv[i+1]);
                break;
            case 'S':
                arg_stream >> switch_tol;
                // switch_tol = strtod(argv[i+1],0);
                break;
            case 'i':
                arg_stream >> increment;
                // increment = strtod(argv[i + 1], 0);
                break;
            case 'I':
                Gillespy::map_variable_populations(arg_stream);
                break;
            case 'p':
                Gillespy::map_variable_parameters(arg_stream);
                break;
            case 'T':
                arg_stream >> trajectories;
                // printf("trajectories: %d\n", trajectories);
                // trajectories = atoi(argv[i+1]);
                break;
            case 0:
                arg_stream >> tau_tol;
                // tau_tol = strtod(argv[i + 1], 0);
                break;
            default:
                std::cerr << "error in option matching" << std::endl;
                break;
        }        
    }

}

ArgParser::~ArgParser()
{
    // delete args;
}
// int main(int argc, char *argv[])
// {
//     ArgParser parser = ArgParser(argc, argv);
//     std::cout << parser.get("timesteps");
//     return 0;
// }
