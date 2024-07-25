#ifndef POLYGONMATCHING_COMMAND_LINE_PARSER_H
#define POLYGONMATCHING_COMMAND_LINE_PARSER_H

#include <string>
#include <optional>
#include <stdexcept>

//struct to store set command line options
struct CommandLineOptions {
    std::string dataset_name;
    double lambda = -1.0; // set default to -1.0 to recognize non-set lambda
    int num_threads = 1; //set default to 1 thread
    int time_limit = 1000; //set default to 1000s for the time limit per component
    int size_limit = -1; //size limit for sets to me matched, deactivate on default
};

void print_help();
CommandLineOptions parse_command_line(int argc, char* argv[]);

#endif //POLYGONMATCHING_COMMAND_LINE_PARSER_H
