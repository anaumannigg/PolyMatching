#include "../include/command_line_parser.h"

#include <iostream>
#include <sstream>

void print_help() {
    std::cout << "Usage: polygonmatching [-d dataset_name] [-l lambda]\n"
              << "Options:\n"
              << "  -d dataset_name      Provide the name of the dataset. There should be two Shapefiles\n"
              << "                       of SINGLE(!) Polygons with appendices _osm and _alkis.\n"
              << "                       The code expects it to be located in 'input/dataset_name'.\n"
              << "  -l lambda            Specify the parameter for lambda.\n"
              << "  -t threads           (optional, default 1) Specify the number of threads that will be used for computations on connected components.\n"
              << "  -k time_limit        (optional, default is 1000) Specify the time limit in [s] for the computation per connected component.\n"
              << "  -s size_limit        (optional, deactivated on default) Set a limit on the set size on cumulative vertex computation.\n"
              << "                       Note: m:n-matches might still be bigger due to precomputations, this only affects the graph traversal.\n"
              << "  -h, --help           Display this help message.\n";
}

CommandLineOptions parse_command_line(int argc, char* argv[]) {
    CommandLineOptions options;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            print_help();
            exit(0);
        } else if (arg == "-d" && i + 1 < argc) {
            options.dataset_name = argv[++i];
        } else if (arg == "-l" && i + 1 < argc) {
            options.lambda = atof(argv[++i]);
        } else if (arg == "-t" & i + 1 < argc) {
            options.num_threads = atoi(argv[++i]);
        } else if (arg == "-k" & i + 1 < argc) {
            options.time_limit = atoi(argv[++i]);
        } else if (arg == "-s" & i + 1 < argc) {
            options.size_limit = atoi(argv[++i]);
        } else {
            throw std::invalid_argument("Unknown option: " + arg);
        }
    }
    return options;
}