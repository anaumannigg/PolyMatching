#include "../include/cgal_includes.h"
#include "../include/shapefile_io_operations.h"
#include "../include/localization.h"
#include "../include/graph_computations.h"
#include "../include/linear_program.h"
#include "../include/command_line_parser.h"
#include "../include/binary_io.h"

#include <exception>
#include <filesystem>
namespace fs = std::filesystem;

//a thread for printing out the status of computation
void statusTHREAD(std::atomic<int>& processed_counter, int num_sets) {
    while (true) {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        int processed = processed_counter.load();
        // Move the cursor up one line and clear it
        std::cout << "\x1b[1A\x1b[2K";
        std::cout << "Processed: " << processed << " / " << num_sets << " sets." << std::endl;
        if (processed >= num_sets) break; // Exit when all sets are processed
    }
}

void setDecompositionTHREAD(const std::vector<Polygon_wh>& osm_polys, const std::vector<Polygon_wh>& atkis_polys, const std::vector<Polygon_wh>& merged_polys,
                            std::pair<int,int> interval, int thread_id, const Localization& osm_rtree, const Localization& atkis_rtree,
                            std::vector<std::vector<int>>& osm_poly_id_lookup, std::vector<std::vector<int>>& atkis_poly_id_lookup, std::vector<int>& set_sizes,
                            std::string folderMERGED, std::string data_name, std::atomic<int>& processed_counter) {

    //specifiy path for binary file
    //since the number of decomposition threads is equal to the number of worker threads later,
    //each decomposition thread creates exactly one binary file for a worker thread to read later on
    std::string binary_path = folderMERGED + "data_thread" + to_string(thread_id);

    //ensure file does not exist yet to not append to already existing file
    std::ifstream file(binary_path);
    if (file.good()) {
        file.close();
        std::remove(binary_path.c_str());
    }

    for (int p = interval.first; p <= interval.second; p++) {
        const Polygon_wh m_poly = merged_polys[p];
        //collect all intersecting polygons in both maps
        std::vector<int> osm_group_indices, atkis_group_indices;
        osm_rtree.get_neighbors(m_poly, &osm_group_indices);
        atkis_rtree.get_neighbors(m_poly, &atkis_group_indices);

        //remember indices for lookup
        osm_poly_id_lookup[p] = osm_group_indices;
        atkis_poly_id_lookup[p]= atkis_group_indices;

        std::vector<Polygon_wh> osm_group, atkis_group;
        for (const auto& i : osm_group_indices) osm_group.push_back(osm_polys[i]);
        for (const auto& i : atkis_group_indices) atkis_group.push_back(atkis_polys[i]);

        //store size of set
        set_sizes[p] = osm_group.size() + atkis_group.size();

        //write to binary_file
        writePolysToBinaryFile(osm_group,atkis_group,binary_path,p);

        processed_counter++;
    }
}

//run n : m map matching algorithm on the maps provided by fileOSM and fileATKIS for each of the provided epsilon values
//returns Solution to append to global Solution
void solveConnectedSetTHREAD(std::string data_name, double lambda, int time_limit, int size_limit, std::pair<int,int> set_intervals, int thread_id,
                             std::vector<bool>& completed_exploration, std::vector<std::pair<int,int>>& set_sizes,
                             std::vector<Solution>& sols, std::vector<double>& execution_times, std::atomic<int>& processed_counter) {

    //get path of binary file for this thread
    std::string binary_path = "../input/" + data_name + "/merged/data_thread" + std::to_string(thread_id);

    //create binary reader object, which will parse the bin file line by line
    BinaryPolygonFileReader bin_reader(binary_path);

    //init set id field
    size_t set_id = 0;
    std::vector<Polygon_wh> osm_polys,atkis_polys;

    //this thread should handle the sets given in the range of indices provided
    while(bin_reader.readNextSet(set_id, osm_polys,atkis_polys)) {
        //measure time for anaylsis output
        std::chrono::steady_clock::time_point set_begin = std::chrono::steady_clock::now();

        //collect and output all matched unmatched polygons from osm and atkis polys
        std::vector <Polygon_wh> osm_unconsidered, osm_considered;
        std::vector <Polygon_wh> atkis_unconsidered, atkis_considered;

        /*
        //DATA CLEANING: make sure no polygons of one set intersect each other (this might happen due to artifacts
        //in data sets such as underground levels)
        for (std::vector <Polygon_wh> &poly_set: {std::ref(osm_polys), std::ref(atkis_polys)}) {
            for (int i = 0; i < poly_set.size(); i++) {
                for (int j = i + 1; j < poly_set.size(); j++) {
                    if (CGAL::do_intersect(poly_set[i], poly_set[j])) {
                        //get intersection area
                        std::list <Polygon_wh> intersection;
                        CGAL::intersection(poly_set[i], poly_set[j], std::back_inserter(intersection));

                        double inter_area = 0.0;
                        for (const auto &i_poly: intersection) {
                            if (i_poly.outer_boundary().is_simple()) {
                                inter_area += to_double(i_poly.outer_boundary().area());
                            }
                        }
                        //if intersection area is not negligible, delete the bigger polygon
                        if (inter_area > 0.1) {
                            double i_area = to_double(poly_set[i].outer_boundary().area());
                            double j_area = to_double(poly_set[j].outer_boundary().area());
                            if (i_area > j_area) {
                                poly_set.erase(std::next(poly_set.begin(), i));
                                i--;
                                j = poly_set.size();
                            } else {
                                poly_set.erase(std::next(poly_set.begin(), j));
                                j--;
                            }
                        }
                    }
                }
            }

        }
        //END DATA CLEANING*/


        //init r-trees of osm and atkis datasets for spatial queries
        Localization osm_rtree(osm_polys);
        Localization atkis_rtree(atkis_polys);

        //precompute all areas within the arrangement of lines formed by the polygons in order to save computation time
        MapOverlay mo = MapOverlay(osm_polys, atkis_polys);
        mo.assignFaces(osm_polys, atkis_polys, osm_rtree, atkis_rtree);

        //initialize sol_local storage
        Solution sol_local(osm_polys.size(), atkis_polys.size());

        //initialize vector of graphs for sol_local on threads
        std::vector <CandidateGraph> g_vec;

        //init lists to remember, which polygons were already considered
        std::vector<bool> osm_poly_visited;
        osm_poly_visited.resize(osm_polys.size(), false);
        std::vector<bool> atkis_poly_visited;
        atkis_poly_visited.resize(atkis_polys.size(), false);

        int max_graph_size = 0;
        int number_of_cons_vertices = 0;

        //iterate over OSM polygons and find neighbors
        int osm_index = 0;
        for (const auto &osm_p: osm_polys) {

            //if polygon is already visited_continue
            if (osm_poly_visited[osm_index]) {
                osm_index++;
                continue;
            } else osm_poly_visited[osm_index] = true;

            //build multilayer graph of the connected component including the osm polygon
            CandidateGraph cg;

            //add first node to graph
            std::vector<int> v_polys = {osm_index};
            int vertex_zero = cg.add_vertex(false, v_polys);

            //while new neighbors are found, iteratively collect and add neighbors in both maps
            //map switch variable switches between each loop (0 = OSM, 1 = ATKIS)
            bool map_switch = true;

            //remember neighbors from last loop
            std::vector<int> query_input;
            query_input.push_back(osm_index);

            //begin exploring the dataset and building the intersection graph
            do {
                //vector to store all found neighbors
                std::vector<int> neighbors;
                //vector to store all neighbors, which are fully included in the same opposing polygon and thus should be cumulated already
                std::vector <std::vector<int>> neighbors_grouped;

                //decide, in which map to queue
                Localization *rtree = !map_switch ? &osm_rtree : &atkis_rtree;
                std::vector <Polygon_wh> polys = map_switch ? osm_polys : atkis_polys;

                //queue each polygon in query input
                for (const auto &q: query_input) {
                    //find all neighbors of query input q
                    std::vector<int> neighbors_all;
                    rtree->get_neighbors(polys[q], &neighbors_all);

                    //find all fully included neighbors, which should be introduced as a cumulative vertex right away
                    std::vector<int> neighbors_fully_included;
                    rtree->get_neighbors_fully_included(polys[q], &neighbors_fully_included);

                    //do not consider the fully included neighbors anymore
                    std::sort(neighbors_all.begin(), neighbors_all.end());
                    std::sort(neighbors_fully_included.begin(), neighbors_fully_included.end());

                    //if more than one neighbor is fully included, delete them from the all neighbors set as they should be a cumulated node right away
                    if (neighbors_fully_included.size() > 0) {
                        for (const auto &nfi: neighbors_fully_included)
                            neighbors_all.erase(std::remove(neighbors_all.begin(), neighbors_all.end(), nfi),
                                                neighbors_all.end());
                    }

                    //add to neighbors and cumulative neighbor sets
                    neighbors.insert(neighbors.end(), neighbors_all.begin(), neighbors_all.end());
                    //EXPERIMENT:: also consider single fully included vertices as group
                    if (neighbors_fully_included.size() > 0) {
                        neighbors_grouped.push_back(neighbors_fully_included);
                    }


                }

                //make neighbor vector unique
                std::sort(neighbors.begin(), neighbors.end());
                neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
                std::sort(neighbors_grouped.begin(), neighbors_grouped.end());
                neighbors_grouped.erase(std::unique(neighbors_grouped.begin(), neighbors_grouped.end()),
                                        neighbors_grouped.end());

                //delete all grouped polygons from the single neighbors list, duplicates may happen due to exploration of multiple polygons beforehand
                for (const auto &n_group: neighbors_grouped) {
                    for (const auto &n: n_group)
                        neighbors.erase(std::remove(neighbors.begin(), neighbors.end(), n), neighbors.end());
                }

                //empty query input list for next loop
                query_input.clear();

                //remember, if any new nodes were added to the graph in this iteration
                bool added_new_vertex = false;
                //filter neighbors to only those, which were not found yet, also remember their indices within the corresponding list
                for (const auto &n: neighbors) {
                    std::vector<bool> *search_list = !map_switch ? &osm_poly_visited : &atkis_poly_visited;
                    if (!(*search_list)[n]) {
                        //mark as visited
                        (*search_list)[n] = true;
                        //save vertex info
                        std::vector<int> v_polys = {n};
                        //neighbor polygon has not been considered yet, create a node for it in the graph and set vertex info
                        int new_vertex = cg.add_vertex(map_switch, v_polys);

                        //remember query input for next loop
                        query_input.push_back(n);

                        added_new_vertex = true;
                    }
                }

                //filter grouped neighbors to only those, which were not found yet, also remember their indices within the corresponding list
                for (const auto &n_group: neighbors_grouped) {
                    std::vector<bool> *search_list = !map_switch ? &osm_poly_visited : &atkis_poly_visited;

                    //check if at least one of the polygons in the set is not visited yet
                    bool completely_visited = true;
                    for (const auto &n: n_group) {
                        if (!(*search_list)[n]) {
                            completely_visited = false;
                            break;
                        }
                    }

                    if (!completely_visited) {
                        //mark as visited
                        for (const auto &n: n_group) (*search_list)[n] = true;
                        //save vertex info
                        std::vector<int> v_polys = n_group;
                        //neighbor polygon has not been considered yet, create a node for it in the graph and set vertex info
                        int new_vertex = cg.add_vertex(map_switch, v_polys);

                        //remember query input for next loop
                        for (const auto &n: n_group) {
                            //EXPERIMENT: add another property: only further explore polygons that have another polygon that is at least 50% included in them
                            Localization *rtree_opp = map_switch ? &osm_rtree : &atkis_rtree;
                            std::vector <Polygon_wh> polys_opp = !map_switch ? osm_polys : atkis_polys;
                            std::vector<int> neighbors;
                            rtree_opp->get_neighbors_majorly_included(polys_opp[n], &neighbors, 0.5);
                            //only consider polygon for further exploration if there is a neighbor fulfilling the property
                            if (neighbors.size() > 1) {
                                query_input.push_back(n);
                            }
                        }

                        added_new_vertex = true;
                    }

                }

                //break, once no more neighbors were found
                if (!added_new_vertex) break;

                map_switch = !map_switch;
            } while (true);


            Graph *g = cg.get_graph();
            if (g->vertex_set().size() > max_graph_size) max_graph_size = g->vertex_set().size();
            number_of_cons_vertices += g->vertex_set().size();

            //next step: build cumulated graph with nodes representing sets of polygons
            //if graph is small anyway, do not get into cumulated graph computation
            if (g->vertex_set().size() <= 2) {
                //found isolated vertex, insert it into the sol_local as isolated
                std::vector<int> o;
                o.push_back(osm_index);
                std::vector<int> a;

                if (g->vertex_set().size() == 1) {
                    //do not add a match in this case, as singular polygon matches should be assigned negative indices via "completeMatching" in the end
                    //remember osm node as isolated node for output
                    osm_unconsidered.push_back(osm_polys[(*g)[0].referenced_polys[0]]);
                } else {
                    //graph has size 2, atkis poly should also be pushed back to sol_local
                    //poly index has to be vertex 1, as vertex 0 is always the initial OSM poly
                    for (const auto &atkis_index: (*g)[1].referenced_polys) a.push_back(atkis_index);

                    double q = mo.getIoU(o, a) - lambda;

                    //only form a match, if it leads to a better sol_local quality
                    if (q > 0) sol_local.addMatch(o, a, q);

                }
            } else {
                //graph has size >= 3, remember for thread traversal
                g_vec.push_back(cg);

            }

            osm_index++;
        }

        //only need to compute edge weights, cumulative vertices and solve ILP, if there are graphs with size > 3
        //those are now stored within g_vec
        if (g_vec.size() > 0) {
            //first compute all inclusion properties between polygons
            setInclusionProperties(g_vec, osm_rtree, atkis_rtree, mo, lambda, osm_polys, atkis_polys);

            //next perform the pregrouping, where polygons of the same map, that are always in the same match,
            //get grouped into a single cumulative vertex
            int n_vertices_after_pregrouping = 0;
            for (auto &g: g_vec) {
                pregroupIncludedPolygons(&g, osm_polys, atkis_polys, osm_rtree, atkis_rtree, mo, lambda);
                n_vertices_after_pregrouping += g.num_vertices();
            }
            //remember total vertex count after pregrouping for analysis
            set_sizes[set_id].first = n_vertices_after_pregrouping;


            //compute edges of each graph
            computeEdges(g_vec, osm_rtree, atkis_rtree, mo, lambda,  osm_polys, atkis_polys);


            int n_vertices_after_simple_matches = 0;

            //start precomputation of simple 1:1 matches, which can be excluded from the graph
            for (auto &g: g_vec) {
                precomputeSimpleMatches(&g, mo, lambda, &sol_local);
                n_vertices_after_simple_matches += g.num_vertices();
            }

            //remember total amount of vertices after simple match precomputation for analysis
            set_sizes[set_id].second = n_vertices_after_simple_matches;

            //precomputation might have led to disconnected graphs, which can be deconstructed into connected components again
            //create temporary memory for new g_threads vector
            std::vector <CandidateGraph> g_vec_decomposed;
            for (int g_id = 0; g_id < g_vec.size(); g_id++) {
                Graph g = *g_vec[g_id].get_graph();

                //skip for empty graphs
                if (num_vertices(g) == 0) continue;

                //store component per vertex
                std::vector<int> comps(num_vertices(g));

                //compute connected components
                int num_components = connected_components(g, &comps[0]);

                //if graph can be decomposed into more than one component, it should be split up
                if (num_components > 1) {
                    //store copies of graph
                    std::vector <CandidateGraph> g_comps;
                    for (int c = 0; c < num_components; c++) g_comps.push_back(g_vec[g_id].copy());

                    //loop down to not change vertex numbers
                    for (int v = comps.size() - 1; v >= 0; v--) {
                        //the graphs in g_comps should represent the connected components, delete vertex from every graph except the one representing the component it is included in
                        for (int c = 0; c < num_components; c++) {
                            if (comps[v] != c) {
                                g_comps[c].delete_vertex(v);
                            }

                        }
                    }

                    //add new graphs to vector
                    for (int c = 0; c < num_components; c++) {
                        if (g_comps[c].num_vertices() > 1) g_vec_decomposed.push_back(g_comps[c]);
                    }


                } else {
                    //graph stays the same, insert to new list while filtering out empty graphs or those with only one vertex that might occur due to prematching simple matches
                    if (g_vec[g_id].num_vertices() > 1) g_vec_decomposed.push_back(g_vec[g_id]);
                }
            }

            //reset g_threads pointer
            g_vec.clear();
            g_vec = g_vec_decomposed;

            //start taking time for eventually killing the threads if they take too long
            std::chrono::steady_clock::time_point components_begin = std::chrono::steady_clock::now();

            //start threads computing the candidate graph
            for (auto& g : g_vec) {
                completeCandidateGraph(&g, mo, lambda, std::chrono::steady_clock::now(), time_limit, size_limit);
            }

            //after component computation terminated, check if time exceeded the threshold, if so mark as incomplete
            std::chrono::steady_clock::time_point components_end = std::chrono::steady_clock::now();
            if (components_end - components_begin >= std::chrono::seconds(time_limit)) {
                std::cout << "\x1b[1A\x1b[2K";
                cout << "Time limit reached on set " << set_id << ".\n\n";
                completed_exploration[set_id] = false;
            } else completed_exploration[set_id] = true;

            //for each graph set up ILP to solve cumulated graph problem optimally (GUROBI is already multi-threaded)
            for (auto graph: g_vec) {
                LinearProgram::solveILP(*graph.get_graph(), osm_polys.size(), atkis_polys.size(), &sol_local);
            }

        }


        //complete the matching in the sol_local via assigning unmatched polygons negative indices
        sol_local.completeMatching();


        sols[set_id] = sol_local;

        std::chrono::steady_clock::time_point set_end = std::chrono::steady_clock::now();
        execution_times[set_id] = std::chrono::duration_cast<std::chrono::milliseconds>(set_end - set_begin).count();

        processed_counter++;
    }

}

//takes two maps and a third map that contains the merged version of both maps with a union of all polygons with a common boundary
//decomposes the maps w.r.t. the merged map into separated subsets and saves them in a subfolder 'tmp' of the current folder
void run_NtoM_with_PreDecomposition(const CommandLineOptions& options) {
    //extract options
    std::string data_name = options.dataset_name;
    double lambda = options.lambda;
    int num_threads = options.num_threads;
    int time_limit = options.time_limit;
    int size_limit = options.size_limit;

    //timekeeping
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    //retreive paths
    std::string fileOSM = "../input/" + data_name + "/" + data_name + "_osm";
    std::string fileATKIS = "../input/" + data_name + "/" + data_name + "_atkis";


    std::string folderMERGED = "../input/" + data_name + "/merged/";
    //create folder of it does not exist yet
    if (!fs::exists(folderMERGED)) fs::create_directory(folderMERGED);

    //read OSM polygons from shapefile
    std::vector <Polygon_wh> osm_polys, atkis_polys, merged_polys;
    ReadShapeFile(fileOSM, &osm_polys);
    ReadShapeFile(fileATKIS, &atkis_polys);

    cout << "read " << osm_polys.size() << " polygons from osm file." << endl;
    cout << "read " << atkis_polys.size() << " polygons from atkis file." << endl;

    //if a file with the merged information is already provided, load it, else compute the union and save it for
    //later use
    std::string fileMERGED = "../input/" + data_name + "/" + data_name + "_merged";
    if (fs::exists(fileMERGED + ".shp")){
        ReadShapeFile(fileMERGED, &merged_polys);
        cout << "read " << merged_polys.size() << " polygons from merged file." << endl;

    }
    else {
        merged_polys = merge(osm_polys,atkis_polys);
        writeToShapeFile(merged_polys, fileMERGED);
        cout << "computed " << merged_polys.size() << " components via merging." << endl;
    }

    //check if any polygon within the sets is non-simple, as this will lead to issues in the code
    //warn the user if a non-simple polygon is found
    for(const auto& polyset: {osm_polys,atkis_polys,merged_polys}) {
        for(const auto& p : polyset)  {
            if(!p.outer_boundary().is_simple()) {
                cout << "WARNING: Found non-simple polygon in the input, that will lead to errors in the code: ";
                cout << "POLYGON((";
                for(const auto& v : p.outer_boundary().vertices()) cout << std::fixed << v << ", ";
                cout<<"))\n";
            }
        }
    }

	//init r-trees of osm and atkis datasets in order for quicker neighbor finding
	Localization osm_rtree(osm_polys);
	Localization atkis_rtree(atkis_polys);

	//init original vertex ID lookup, as IDs change on grouping
    int num_components = merged_polys.size();
	std::vector<std::vector<int>> osm_poly_id_lookup(num_components), atkis_poly_id_lookup(num_components);

    //set up intervals for the computation on the connected sets using multiple threads
    //the intervals will be used for the decomposition as well as the solving
    num_threads = std::min(num_threads, num_components);
    int sets_per_thread = (int)(num_components / num_threads);

    std::vector<std::pair<int,int>> thread_set_intervals;
    for(int t=0; t<num_threads-1;t++) {
        thread_set_intervals.emplace_back(t * sets_per_thread, (t + 1) * sets_per_thread - 1);
    }
    thread_set_intervals.emplace_back((num_threads-1) * sets_per_thread, num_components - 1);


    //loop over all merged polygons and create sets (<-> connected components in the intersection graph)
    std::vector<int> set_sizes(num_components);

    cout << "Starting Decomposition into connected sets.\n\n";

    //start status thread
    std::atomic<int> processed_counter(0);
    std::thread status_thread(statusTHREAD, std::ref(processed_counter),num_components);

    //start worker threads
    std::vector<std::thread> decomp_threads(num_threads);
    for(int t=0; t < num_threads; t++) {
        decomp_threads[t] = std::thread(setDecompositionTHREAD,
                                     std::ref(osm_polys), std::ref(atkis_polys),std::ref(merged_polys),
                                     thread_set_intervals[t], t, std::ref(osm_rtree), std::ref(atkis_rtree),
                                     std::ref(osm_poly_id_lookup), std::ref(atkis_poly_id_lookup), std::ref(set_sizes),
                                     folderMERGED, data_name, std::ref(processed_counter));
    }

    //join threads
    for(int t=0; t < num_threads; t++) {
        decomp_threads[t].join();
    }
    //ensure the status thread finished
    processed_counter = num_components;
    status_thread.join();

	cout << "Completed inital decomposition into connected subsets." << endl;

	//initialize global solution
	Solution sol(osm_polys.size(),atkis_polys.size());

	//remember for each instance, if the exploration could be performed within the time limits
	std::vector<bool> completed_exploration; completed_exploration.resize(num_components, true);
	std::vector<double> execution_times; execution_times.resize(num_components, -1.0);




    std::vector<Solution> thread_solutions(num_components);
    std::vector<bool> thread_explored(num_components);
    std::vector<std::pair<int,int>> set_sizes_after_precomputations(num_components);

    cout << "Starting to solve sets.\n\n";

    //start status thread
    processed_counter = 0;
    status_thread = std::thread(statusTHREAD, std::ref(processed_counter), num_components);

    //start worker threads
    std::vector<std::thread> set_threads(num_threads);
    for(int t=0; t < num_threads; t++) {
        set_threads[t] = std::thread(solveConnectedSetTHREAD,
                                     data_name, lambda, time_limit, size_limit, thread_set_intervals[t], t,
                                     std::ref(thread_explored), std::ref(set_sizes_after_precomputations),
                                     std::ref(thread_solutions), std::ref(execution_times), std::ref(processed_counter));
    }

    //join threads
    for(int t=0; t < num_threads; t++) {
        set_threads[t].join();
    }
    // Ensure the status thread finishes
    processed_counter = num_components;
    status_thread.join();

    //insert all part solutions
    for(int s=0; s < num_components; s++) {
        sol.insert(thread_solutions[s], osm_poly_id_lookup[s], atkis_poly_id_lookup[s]);
    }


	//complete global matching (setting negative match IDs for unmatched polygons)
	sol.completeMatching();

	//make sure the target folder exists
	std::string target_folder = "../export/" + data_name + "/lambda" + to_string((int)(lambda * 100));
	if(!fs::exists(target_folder)) fs::create_directories(target_folder);

	//export the global solution
	std::string osm_output_path = "../export/" + data_name + "/lambda" + to_string((int)(lambda * 100)) + "/" + data_name + "_osm_matched_NtoM";
	std::string atkis_output_path = "../export/" + data_name + "/lambda" + to_string((int)(lambda * 100)) + "/" + data_name + "_atkis_matched_NtoM";
	std::string csv_output_path = "../export/" + data_name + "/lambda" + to_string((int)(lambda * 100)) + "/" + data_name + "_data_NtoM";
	std::string analysis_output_path = "../export/" + data_name + "/lambda" + to_string((int)(lambda * 100)) + "/" + data_name + "_analysis";

    //write labeled polygons to shapefiles
	writeToShapeFile(osm_polys, sol.getMatchIndices(0), sol.getMatchWeights(0), osm_output_path);
	writeToShapeFile(atkis_polys, sol.getMatchIndices(1), sol.getMatchWeights(1), atkis_output_path);
	//write analysis data to csvs
	writeToCSV(sol, csv_output_path);
	writeToCSV(set_sizes, set_sizes_after_precomputations, execution_times, completed_exploration, analysis_output_path);


	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	cout << "Computed n:m matching and exported results. (" << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << " [s])" << endl;
	cout << "Objective = " << sol.getTargetValue() << endl;

}

int main(int argc, char* argv[]) {

    try {
        CommandLineOptions options = parse_command_line(argc,argv);

        if (options.dataset_name.empty() ||options.lambda == -1.0) {
            std::cerr << "Error: Both -d and -l options are required." << endl;
            print_help();
            return 1;
        }

        run_NtoM_with_PreDecomposition(options);

        cout << "Completed matching " << options.dataset_name << " for lambda " << std::to_string(options.lambda) << ".\n";

    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << endl;
        print_help();
        return 1;
    }


	return 0;
}
