#include <map>
#include <string>
#include <iostream>
#include <getopt.h>

class ArgParser
{
private:
    std::map<char*, int> args; 
    std::map<char*, char*> parameters;
    std::vector<int> init_pop;
    // int init_pop[];
public:
    int get(char* arg) {
        return args[arg];
    };
    ArgParser(const int argc, char *argv[]);
    ~ArgParser();
};

ArgParser::ArgParser(const int argc, char* argv[])
{
    args = std::map<char*, int>();
    int opt;
    while (1) {
        static struct option long_options[] = {
            {"timesteps", required_argument, 0, 't'},
            {"end", required_argument, 0, 'e'},
            {"seed", required_argument, 0, 's'},
            {"trajectories", required_argument, 0, 'T'},
            {"increment", required_argument, 0, 'i'},
            {"tau_tol", required_argument, 0, 0},
            {"switch_tol", required_argument, 0, 'S'},
            {"initpop", required_argument, 0, 'I'},
            {"parameters", required_argument, 0, 'P'},
            {0,0,0,0}
        };
        opt = getopt_long(argc, argv, "t:e:s:T:i:0:S:I:P:",long_options, NULL);
        // std::cout << opt;
        if (opt == -1)
            break;
        
        switch (opt)
        {
        case 0:
            args.insert(std::pair<char *, int>("tau_tol", atoi(optarg)));
            break;
        case 't':
            args.insert(std::pair<char *, int>("timesteps", atoi(optarg)));
            break;
        case 'e':
            args.insert(std::pair<char *, int>("end", atoi(optarg)));
            break;
        case 's':
            args.insert(std::pair<char *, int>("seed", atoi(optarg)));
            break;
        case 'T':
            args.insert(std::pair<char *, int>("trajectories", atoi(optarg)));
            break;
        case 'i':
            args.insert(std::pair<char *, int>("increment", atoi(optarg)));
            break;
        case 'S':
            break;
        case 'I':
            // args.insert(std::pair<char *, int>("initpop", atoi(optarg)));
            break;
        case 'P':
            // args.insert(std::pair<char *, int>("parameters", atoi(optarg)));
            break;
        
        default:
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
