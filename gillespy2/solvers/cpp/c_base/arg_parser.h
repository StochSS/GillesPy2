#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <getopt.h>
#include <string.h>

class ArgParser
{
private:
    

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
    int trajectories;
    int timesteps;
    double end;
    int seed;
    double increment;
    double switch_tol;
    double tau_tol;
    std::vector<int> init_pop;
    std::vector<double> parameters;
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

ArgParser::ArgParser(int argc, char* argv[]):
    
{
    // long_options 
    // args = std::map<char*, int>();
    int opt;
    // options arg;
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
                end = atoi(argv[i+1]);
                break;
            case 's':
                args["end"] = atoi(argv[i+1]);
                break;
            case 'S':
                break;
            case 'i':
                break;
            case 'I':
                break;
            case 'p':
                break;
            case 'T':
                break;
            case 0:
                break;
            default:
                break;        }
    }
    // while (1) {
    //     static struct option long_options[] = {
    //         {"timesteps", required_argument, 0, 't'},
    //         {"end", required_argument, 0, 'e'},
    //         {"seed", required_argument, 0, 's'},
    //         {"trajectories", required_argument, 0, 'T'},
    //         {"increment", required_argument, 0, 'i'},
    //         {"tau_tol", required_argument, 0, 0},
    //         {"switch_tol", required_argument, 0, 'S'},
    //         {"initpop", required_argument, 0, 'I'},
    //         {"parameters", required_argument, 0, 'P'},
    //         {0,0,0,0}
    //     };
    //     opt = getopt_long(argc, argv, "t:e:s:T:i:0:S:I:P:",long_options, NULL);
    //     // std::cout << opt;
    //     if (opt == -1)
    //         break;
        
    //     switch (opt)
    //     {
    //     case 0:
    //         args.insert(std::pair<char *, int>("tau_tol", atoi(optarg)));
    //         break;
    //     case 't':
    //         args.insert(std::pair<char *, int>("timesteps", atoi(optarg)));
    //         break;
    //     case 'e':
    //         args.insert(std::pair<char *, int>("end", atoi(optarg)));
    //         break;
    //     case 's':
    //         args.insert(std::pair<char *, int>("seed", atoi(optarg)));
    //         break;
    //     case 'T':
    //         args.insert(std::pair<char *, int>("trajectories", atoi(optarg)));
    //         break;
    //     case 'i':
    //         args.insert(std::pair<char *, int>("increment", atoi(optarg)));
    //         break;
    //     case 'S':
    //         break;
    //     case 'I':
    //         // args.insert(std::pair<char *, int>("initpop", atoi(optarg)));
    //         break;
    //     case 'P':
    //         // args.insert(std::pair<char *, int>("parameters", atoi(optarg)));
    //         break;
        
    //     default:
    //         break;
    //     }
    // }

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
