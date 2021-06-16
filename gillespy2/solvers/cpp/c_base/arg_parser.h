#include <map>
#include <string>
#include <vector>
#include <iostream>
// #include <getopt.h>
#include <string.h>
#include <stdlib.h>

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
[--tau_tol] <double> \n\
";
    char match_arg(char* opt)
    {
        if (!strcmp(opt, "timesteps"))
            return 't';
        if (!strcmp(opt, "end"))
            return 'e';
        if (!strcmp(opt, "seed"))
            return 's';
        if (!strcmp(opt, "switch_tol"))
            return 'S';
        if (!strcmp(opt, "increment")) //double
            return 'i';
        if (!strcmp(opt, "init_pop"))
            return 'I';
        if (!strcmp(opt, "parameters"))
            return 'p';
        if (!strcmp(opt, "trajectories"))
            return 'T';
        if (!strcmp(opt, "tau_tol"))
            return 0;      
    };
public:
    int trajectories = 0;
    int timesteps = 0;
    double end = 0.0;
    int seed = -1;
    double increment = 0.0;
    double switch_tol = 0.0;
    double tau_tol = 0.0;
    std::vector<int> init_pop = std::vector<int>();
    std::vector<double> parameters = std::vector<double>();
    // int get(char* arg) {
    //     if (!strcmp(arg, "timesteps"))
    //         return 't';
    //     if (!strcmp(arg, "end"))
    //         return 'e';
    //     if (!strcmp(arg, "seed"))
    //         return 's';
    //     if (!strcmp(arg, "switch_tol"))
    //         return 'S';
    //     if (!strcmp(arg, "increment")) //double
    //         return 'i';
    //     if (!strcmp(arg, "init_pop"))
    //         return 'I';
    //     if (!strcmp(arg, "parameters"))
    //         return 'p';
    //     if (!strcmp(arg, "trajectories"))
    //         return 'T';
    //     if (!strcmp(arg, "tau_tol"))
    //         return 0;
    //     return args[arg];
    // };
    ArgParser(const int argc, char *argv[]);
    ~ArgParser();
};

ArgParser::ArgParser(int argc, char* argv[])    
{
    int opt;
    for (int i = 1; i < argc; ++i){
        if (argv[i][0] == '-' && argv[i][1] != '-') {
            opt = argv[i][1];
        }
        if (argv[i][0] == '-' && argv[i][1] == '-') {
            opt = match_arg((char*)argv[i][2]);
        }
        switch(opt) {
            case 't':
                timesteps = atoi(argv[i+1]);
                break;
            case 'e':
                end = strtod(argv[i+1],0);
                break;
            case 's':
                seed = atoi(argv[i+1]);
                break;
            case 'S':
                switch_tol = strtod(argv[i+1],0);
                break;
            case 'i':
                increment = strtod(argv[i + 1], 0);
                break;
            case 'I':
                for (int j = i+1; argv[j][0]!= '-'; ++j){
                    init_pop.push_back(atoi(argv[j]));
                }
                break;
            case 'p':
                for (int j = i + 1; argv[j][0] != '-'; ++j)
                {
                    parameters.push_back(strtod(argv[j], 0));
                }
                break;
            case 'T':
                trajectories = atoi(argv[i+1]);
                break;
            case 0:
                tau_tol = strtod(argv[i + 1], 0);
                break;
            default:
                std::cerr << usage;
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
