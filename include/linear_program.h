#ifndef _linear_program_included_
#define _linear_program_included_

#include "gurobi_c++.h"
#include "cgal_includes.h"
#include "graph_computations.h"
#include "solution.h"

using namespace std;

class LinearProgram {
	

public:
    //sets up and solves ILP based on a bipartite graph G.
    //maximizes sum of chosen edges in constrained bipartite matching problem
    //the result is written into the provided solution object
	static void solveILP(const Graph& g, int num_osm_polys, int num_atkis_polys, Solution* solution);
};

#endif