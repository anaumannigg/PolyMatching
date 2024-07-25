#include "../include/linear_program.h"

//creates and solves an integer linear program w.r.t the nodes and edge weights within the given Graph g (which is the 2nd layer of the Multilayer-Graph in this case)
void LinearProgram::solveILP(const Graph& g, int num_osm_polys, int num_atkis_polys, Solution* solution) {



    try {
	//gurobi check
	GRBEnv env = GRBEnv(true);

    //deactivate console logging
    env.set("LogToConsole", "0");

    //env.set("LogFile", "mip1.log");
    env.start();

    // Create an empty model
    GRBModel model = GRBModel(env);

    // Create variables per edge in G (variable should be 1 when chosen and 0 when not chosen
    std::vector<GRBVar> var;
    std::vector<double> weight;

    //remember adjacent vertices per edge
    std::vector<std::vector<GRBVar>> var_adj_osm_v(g.vertex_set().size()), var_adj_atkis_v(g.vertex_set().size());
    //remember order of edges
    std::vector<int> sources,targets;


    Graph::vertex_iterator v, vend;
    for (boost::tie(v, vend) = vertices(g); v != vend; ++v) {
        //to not count edges twice, only insert new variables once per outgoing edge of the left side of the graph
        if (!g[*v].referenced_map) {
            for (const auto& e : g.out_edge_list(*v)) {
                GRBVar edge_var = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "");
                var.push_back(edge_var);
                double w = get(boost::edge_weight, e);

                weight.push_back(w);

                //remember edge for atkis and osm polys
                var_adj_osm_v[*v].push_back(edge_var);
                var_adj_atkis_v[e.get_target()].push_back(edge_var);

                //remember edge for solution retreival
                sources.push_back(*v);

                targets.push_back(e.get_target());
            }
        }
    }

    //build target function, which is the sum over all weighted edges
    GRBLinExpr target_func;
    target_func.addTerms(&weight[0], &var[0], var.size());

    //set the objective to maximizing the target function
    model.setObjective(target_func, GRB_MAXIMIZE);

    //add legality constraints
    for (bool map_switch : {false, true}) {
        int num_polys = !map_switch ? num_osm_polys : num_atkis_polys;
        for (int i = 0; i < num_polys; i++) {
            //collect all vertices, which refer to polygon i
            std::vector<int> referring_vertices;
            for (boost::tie(v, vend) = vertices(g); v != vend; ++v) {
                if (g[*v].referenced_map == map_switch && std::find(g[*v].referenced_polys.begin(), g[*v].referenced_polys.end(), i) != g[*v].referenced_polys.end()) {
                    //map fits and i is one of the referenced polygons, remember vertex
                    referring_vertices.push_back(*v);
                }
            }

            //if there are referring vertices in the current map, add expression
            if (referring_vertices.size() > 0) {
                //sum of variables referring to the respective nodes should not be greater than 1, s.t. every polygon gets picked at most once
                GRBLinExpr constr_expr;
                std::vector<std::vector<GRBVar>> adjacent_edges = !map_switch ? var_adj_osm_v : var_adj_atkis_v;
                for (const auto& rv : referring_vertices) {
                    double coeff = 1;
                    for (const auto& adj_edge : adjacent_edges[rv]) {
                        constr_expr.addTerms(&coeff, &adj_edge, 1);
                    }
                }
                //add constraint to model, that the polygon may not be chosen twice
                model.addConstr(constr_expr, GRB_LESS_EQUAL, 1.0, "");
            }

        }
    }

    // Optimize model
    model.optimize();

    int matches_added = 0;
    //retreive solution: for each selected edge, put the referred polygons of source and target into the solution
    for (int i = 0; i < var.size(); i++) {
        if (var[i].get(GRB_DoubleAttr_X) == 1.0) {
            //edge is selected, put pair of source and target into solution including the weight of the connecting edge
            double weight = 0.0; bool weight_found = false;
            for (const auto& e : g.out_edge_list(sources[i])) {
                if ((int)e.get_target() == targets[i]) {
                    weight = get(boost::edge_weight, e);
                    weight_found = true;
                    break;
                }
            }
            assert(weight_found && "did not find weight when iterating through edges!");

            solution->addMatch(g[sources[i]].referenced_polys, g[targets[i]].referenced_polys, weight);

            matches_added++;
        }
    }

    }
    catch (GRBException e) {
        cout << "Gurobi crashed with error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Exception during optimization" << endl;
    }

}