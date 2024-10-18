#include "../include/shapefile_io_operations.h"

using namespace boost::geometry;

//reads shapefile and writes polygon information into argument vector
//stability variables:
//min_poly_point_distance states the minimal distance at which two consecutive points of read polygons will be added, otherwise they will be ignored
//max_poy_angle states the maximum angle at which a new point of a polygon will be added, otherwise they will be ignored
void ReadShapeFile(const std::string& filename, std::vector<Polygon_wh>* polys) {
    try {
        SHPHandle handle = SHPOpen(filename.c_str(), "rb");
        if (!handle) {
            throw std::string("File " + filename + " not found.");
        }

        int nShapeType, nEntities;
        double adfMinBound[4], adfMaxBound[4];
        SHPGetInfo(handle, &nEntities, &nShapeType, adfMinBound, adfMaxBound);

        for (int i = 0; i < nEntities; i++) {
            SHPObject* psShape = SHPReadObject(handle, i);

            //Read Polys with potential holes
            if (psShape->nSHPType == SHPT_POLYGON && psShape->nParts >= 1) {
                //create a polygon for the outer boundary
                Polygon outer_boundary;
                std::vector<Point> already_added;

                //get start and end vertex of outer boundary
                int startV = psShape->panPartStart[0];
                int endV = psShape->nParts == 1 ? psShape->nVertices : psShape->panPartStart[1];

                double* x = psShape->padfX;
                double* y = psShape->padfY;

                //add all points of outer boundary
                for (int v = startV; v < endV; ++v) {
                    Point new_p = Point(x[v],y[v]);
                    if(std::find(already_added.begin(),already_added.end(),new_p) == already_added.end()) {
                        already_added.push_back(new_p);
                        outer_boundary.push_back(Point(x[v], y[v]));
                    }
                }

                /*
                cout << "POLYGON((";
                for(const auto& p : outer_boundary) cout << std::fixed << std::setprecision(3) << p << ", ";
                cout <<"))" << endl;*/

                // Ensure the outer boundary is counterclockwise oriented
                if (!outer_boundary.is_counterclockwise_oriented()) {
                    outer_boundary.reverse_orientation();
                }


                //read all inner holes
                std::vector<Polygon> holes;
                already_added.clear();

                for (int part = 1; part < psShape->nParts; ++part) {
                    Polygon hole;

                    startV = psShape->panPartStart[part];
                    endV = (part == psShape->nParts - 1) ? psShape->nVertices : psShape->panPartStart[part + 1];

                    for (int v = startV; v < endV; ++v) {
                        Point new_p = Point(x[v],y[v]);
                        if(std::find(already_added.begin(),already_added.end(),new_p) == already_added.end()) {
                            already_added.push_back(new_p);
                            hole.push_back(Point(x[v], y[v]));
                        }
                    }

                    // Ensure the hole is clockwise oriented
                    if (hole.is_counterclockwise_oriented()) {
                        hole.reverse_orientation();
                    }

                    holes.push_back(hole);
                    already_added.clear();
                }

                // Create the polygon with holes and add it to the vector
                Polygon_wh poly(outer_boundary, holes.begin(), holes.end());
                polys->push_back(poly);

            }
            SHPDestroyObject(psShape);
        }
        SHPClose(handle);

    }
    catch (const std::string& s) {
        throw s;
    }

}

//writes vector of polygons to shapefile, make sure to add first vertex again at the end as this is the way the read polygon method expects the data to be 
void writeToShapeFile(std::vector<Polygon_wh> polys, std::string path) {

    //create handle
    SHPHandle shapefile = SHPCreate(path.c_str(), SHPT_POLYGON);
    DBFHandle dbfile = DBFCreate(path.c_str());

    //create field for poly ID
    int field_face_id = DBFAddField(dbfile, "ID", FTInteger, 5, 0);

    int f = 0;

    for (const auto& poly : polys) {
        // Count total number of vertices and parts
        int total_vertex_count = 0;
        int part_count = 1 + poly.number_of_holes(); // 1 for outer boundary + holes

        // Get the outer boundary and add its vertex count
        const Polygon& outer_boundary = poly.outer_boundary();
        total_vertex_count += outer_boundary.size();

        // Get hole vertex counts
        for (auto hole_it = poly.holes_begin(); hole_it != poly.holes_end(); ++hole_it) {
            total_vertex_count += hole_it->size();
        }


        //collect vertices
        double* x = new double[total_vertex_count];
        double* y = new double[total_vertex_count];

        // Array to store part start indices
        int* part_start_indices = new int[part_count];
        int vertex_index = 0;

        // Write outer boundary vertices
        part_start_indices[0] = 0; // Outer boundary starts at index 0
        for (auto vit = outer_boundary.vertices_begin(); vit != outer_boundary.vertices_end(); ++vit) {
            x[vertex_index] = to_double(vit->x());
            y[vertex_index] = to_double(vit->y());
            vertex_index++;
        }

        // Write hole vertices and update part start indices
        int part_index = 1;
        for (auto hole_it = poly.holes_begin(); hole_it != poly.holes_end(); ++hole_it) {
            part_start_indices[part_index] = vertex_index;
            //holes need to be stored in clockwise orientation, reverse
            int hole_index = vertex_index + hole_it->size() - 1;
            for (auto vit = hole_it->vertices_begin(); vit != hole_it->vertices_end(); ++vit) {
                x[hole_index] = to_double(vit->x());
                y[hole_index] = to_double(vit->y());
                vertex_index++; hole_index--;
            }
            part_index++;
        }

        //create polygon object
        SHPObject* shape = SHPCreateObject(SHPT_POLYGON, -1, part_count, part_start_indices, nullptr, total_vertex_count, x, y, nullptr, nullptr);

        //write shape into file
        int shape_id = SHPWriteObject(shapefile, -1, shape);

        //set field to face ID
        DBFWriteIntegerAttribute(dbfile, shape_id, field_face_id, f++);

        //memory management
        SHPDestroyObject(shape);
        delete[] x;
        delete[] y;
        delete[] part_start_indices;

    }

    DBFClose(dbfile);
    SHPClose(shapefile);
}

//writes vector of polygons to shapefile, adding the labels stored in the area-vector to the respective polygon
void writeToShapeFile(std::vector<Polygon_wh> polys, std::vector<int> group_index, std::string path) {

    //create handle
    SHPHandle shapefile = SHPCreate(path.c_str(), SHPT_POLYGON);
    DBFHandle dbfile = DBFCreate(path.c_str());


    //create field in database
    int field_match_id = DBFAddField(dbfile, "match", FTInteger, 10, 0);
    int field_poly_num_id = DBFAddField(dbfile, "number", FTInteger, 10, 0);

    int p = 0;

    for (const auto& poly : polys) {
        // Count total number of vertices and parts
        int total_vertex_count = 0;
        int part_count = 1 + poly.number_of_holes(); // 1 for outer boundary + holes

        // Get the outer boundary and add its vertex count
        const Polygon& outer_boundary = poly.outer_boundary();
        total_vertex_count += outer_boundary.size();

        // Get hole vertex counts
        for (auto hole_it = poly.holes_begin(); hole_it != poly.holes_end(); ++hole_it) {
            total_vertex_count += hole_it->size();
        }


        //collect vertices
        double* x = new double[total_vertex_count];
        double* y = new double[total_vertex_count];

        // Array to store part start indices
        int* part_start_indices = new int[part_count];
        int vertex_index = 0;

        // Write outer boundary vertices
        part_start_indices[0] = 0; // Outer boundary starts at index 0
        for (auto vit = outer_boundary.vertices_begin(); vit != outer_boundary.vertices_end(); ++vit) {
            x[vertex_index] = to_double(vit->x());
            y[vertex_index] = to_double(vit->y());
            vertex_index++;
        }

        // Write hole vertices and update part start indices
        int part_index = 1;
        for (auto hole_it = poly.holes_begin(); hole_it != poly.holes_end(); ++hole_it) {
            part_start_indices[part_index] = vertex_index;
            //holes need to be stored in clockwise orientation, reverse
            int hole_index = vertex_index + hole_it->size() - 1;
            for (auto vit = hole_it->vertices_begin(); vit != hole_it->vertices_end(); ++vit) {
                x[hole_index] = to_double(vit->x());
                y[hole_index] = to_double(vit->y());
                vertex_index++; hole_index--;
            }
            part_index++;
        }

        //create polygon object
        SHPObject* shape = SHPCreateObject(SHPT_POLYGON, -1, part_count, part_start_indices, nullptr, total_vertex_count, x, y, nullptr, nullptr);

        //write shape into file
        int shape_id = SHPWriteObject(shapefile, -1, shape);

        //write data into field
        DBFWriteDoubleAttribute(dbfile, shape_id, field_match_id, group_index[p]);
        DBFWriteDoubleAttribute(dbfile, shape_id, field_poly_num_id, p);

        //memory management
        SHPDestroyObject(shape);
        delete[] x;
        delete[] y;
        delete[] part_start_indices;

        p++;

    }

    DBFClose(dbfile);
    SHPClose(shapefile);
}

//writes vector of polygons to shapefile, adding the labels stored in the area- and perimeter-vector to the respective polygon
void writeToShapeFile(std::vector<Polygon_wh> polys, std::vector<int> group_index, std::vector<double> match_weight, std::string path) {
    //create handle
    SHPHandle shapefile = SHPCreate(path.c_str(), SHPT_POLYGON);
    DBFHandle dbfile = DBFCreate(path.c_str());

    //create fields in database
    int field_match_id = DBFAddField(dbfile, "match", FTInteger, 10, 3);
    int field_number_id = DBFAddField(dbfile, "number", FTInteger, 10, 3);
    int field_weight_id = DBFAddField(dbfile, "weight", FTDouble, 10, 3);

    int p = 0;

    for (const auto& poly : polys) {
        // Count total number of vertices and parts
        int total_vertex_count = 0;
        int part_count = 1 + poly.number_of_holes(); // 1 for outer boundary + holes

        // Get the outer boundary and add its vertex count
        const Polygon& outer_boundary = poly.outer_boundary();
        total_vertex_count += outer_boundary.size();

        // Get hole vertex counts
        for (auto hole_it = poly.holes_begin(); hole_it != poly.holes_end(); ++hole_it) {
            total_vertex_count += hole_it->size();
        }


        //collect vertices
        double* x = new double[total_vertex_count];
        double* y = new double[total_vertex_count];

        // Array to store part start indices
        int* part_start_indices = new int[part_count];
        int vertex_index = 0;

        // Write outer boundary vertices
        part_start_indices[0] = 0; // Outer boundary starts at index 0
        for (auto vit = outer_boundary.vertices_begin(); vit != outer_boundary.vertices_end(); ++vit) {
            x[vertex_index] = to_double(vit->x());
            y[vertex_index] = to_double(vit->y());
            vertex_index++;
        }

        // Write hole vertices and update part start indices
        int part_index = 1;
        for (auto hole_it = poly.holes_begin(); hole_it != poly.holes_end(); ++hole_it) {
            part_start_indices[part_index] = vertex_index;
            //holes need to be stored in clockwise orientation, reverse
            int hole_index = vertex_index + hole_it->size() - 1;
            for (auto vit = hole_it->vertices_begin(); vit != hole_it->vertices_end(); ++vit) {
                x[hole_index] = to_double(vit->x());
                y[hole_index] = to_double(vit->y());
                vertex_index++; hole_index--;
            }
            part_index++;
        }

        //create polygon object
        SHPObject* shape = SHPCreateObject(SHPT_POLYGON, -1, part_count, part_start_indices, nullptr, total_vertex_count, x, y, nullptr, nullptr);

        //write shape into file
        int shape_id = SHPWriteObject(shapefile, -1, shape);

        //write data into fields
        DBFWriteIntegerAttribute(dbfile, shape_id, field_match_id, group_index[p]);
        DBFWriteIntegerAttribute(dbfile, shape_id, field_number_id, p);
        DBFWriteDoubleAttribute(dbfile, shape_id, field_weight_id, match_weight[p]);
        p++;

        //memory management
        SHPDestroyObject(shape);
        delete[] x;
        delete[] y;
    }

    DBFClose(dbfile);
    SHPClose(shapefile);
}

void writeToWKT(std::vector<Polygon_wh> polys, std::string path) {

    std::ofstream poly_file;
    poly_file.open(path);

    //insert polygons
    poly_file << "MULTIPOLYGON(";
    for (const auto& poly : polys) {
        poly_file << "((";
        for (const Point v : poly.outer_boundary().vertices()) {
            poly_file << std::fixed << v.x() << " " << v.y() << ",";

        }
        poly_file << ")),\n";
    }
    poly_file << ")" << endl;
    poly_file.close();

}

//writes Segments into File
void writeToWKT(std::vector<Segment> segments, std::string path) {

	std::ofstream poly_file;
	poly_file.open(path);

	//insert polygons
	poly_file << "MULTILINESTRING(";
	for (const auto& s : segments) {
		poly_file << "((";
		for (const Point v :  {s.source(),s.target()}) {
			poly_file << std::fixed << v.x() << " " << v.y() << ",";

		}
		poly_file << ")),\n";
	}
	poly_file << ")" << endl;
	poly_file.close();

}

//writes analysis data to csv file
void writeToCSV(Solution s, std::string path) {
	std::ofstream file;
	file.open(path + ".csv");

	//write header 
	file << "map, poly_id, match_id, match_weight" << endl;

	//insert data per polygon
	for (bool map : {0, 1}) {
		std::vector<double> match_weights = s.getMatchWeights(map);
		std::vector<int> match_indices = s.getMatchIndices(map);
		for (int id = 0; id < match_weights.size(); id++) {
			file << map << "," << id << "," << match_indices[id] << "," << match_weights[id] << endl;
		}
	}
	
	file.close();

}

void writeToCSV(std::vector<int> set_sizes, std::vector<std::pair<int,int>> set_sizes_after_precomp, std::vector<double> execution_times, std::vector<bool> completed_exploration, std::string path) {
	std::ofstream file;
	file.open(path + ".csv");

	//write header 
	file << "set_id,set_size,set_size after pregrouping,set_size after simple matches,execution_time,completed_exploration" << endl;

	//insert data 
	for (int id = 0; id < set_sizes.size(); id++) {
		file << id << "," << set_sizes[id] << "," << set_sizes_after_precomp[id].first << "," << set_sizes_after_precomp[id].second << "," << std::fixed << std::setprecision(3) << execution_times[id] << "," << completed_exploration[id] << endl;
	}

	file.close();

}

//reads out two csv files and puts out indices of matches in fileA, which are different to the matches of fileB
void findDifferingMatches(std::string fileA, std::string fileB) {
	std::ifstream csvA, csvB;
	csvA.open(fileA); csvB.open(fileB);
	//read both files and store as vectors
	std::vector<int> matchesA, matchesB;

	//skip headers
	std::string line; 
	getline(csvA, line);

	while (getline(csvA, line)) {
		//parse line, we want the third line, as this stores the match id
		std::vector<std::string> columns;
		std::stringstream ss(line);
		while(ss.good()) {
			std::string substr;
			getline(ss, substr, ',');
			columns.push_back(substr);
		}


		matchesA.push_back(stoi(columns[2]));
	}

	getline(csvB, line);

	while (getline(csvB, line)) {
		//parse line, we want the third line, as this stores the match id
		std::vector<std::string> columns;
		std::stringstream ss(line);
		while (ss.good()) {
			std::string substr;
			getline(ss, substr, ',');
			columns.push_back(substr);
		}


		matchesB.push_back(stoi(columns[2]));
	}

	std::vector<int> visited_match_ids;

	//remember if any differences have been found
	bool found_difference = false;

	//now the matches per polygon are stores in the matches-vectors
	//iterate over every polygon in fileA 
	for (int i = 0; i < matchesA.size(); i++) {
		int match_id_A = matchesA[i];

		//check if match was already visited
		if (std::find(visited_match_ids.begin(), visited_match_ids.end(), match_id_A) != visited_match_ids.end()) continue;

		//has not been visited yet, add to visited memory
		visited_match_ids.push_back(match_id_A);

		std::vector<int> polys_in_match_A = { i };

		//find every polygon with the same match ID
		for (int j = i+1; j < matchesA.size(); j++) {
			if (matchesA[j] == match_id_A) polys_in_match_A.push_back(j);
		}

		//check the same for the opposing map
		int match_id_B = matchesB[i];
		std::vector<int> polys_in_match_B = { i };
		for (int j = 0; j < matchesB.size(); j++) {
			if (j == i) continue;
			if (matchesB[j] == match_id_B) polys_in_match_B.push_back(j);
		}

		// make sure both vectors are sorted
		std::sort(polys_in_match_A.begin(), polys_in_match_A.end());
		std::sort(polys_in_match_B.begin(), polys_in_match_B.end());

		//if the vectors are not equal, a differing match has been found, output this match
		if (polys_in_match_A != polys_in_match_B) {
			found_difference = true;
			cout << "Found differing match! ID in first matching: " << match_id_A << " | ID in second matching: " << match_id_B << endl;

		}
	}	

	if (!found_difference) cout << "No differences have been found, the two matchings are equal!" << endl;
}

//creates a shapefile including LineStrings that approximate an Arrangement of circular arcs and linearcs
void writeToShapeFile(Arrangement arr, std::string path) {
    int curve_approx_depth = 3;

    //create handle
    SHPHandle shapefile = SHPCreate(path.c_str(), SHPT_ARC);
    DBFHandle dbfile = DBFCreate(path.c_str());

    //create field for arc ID
    int field_face_id = DBFAddField(dbfile, "ID", FTInteger, 5, 0);

    int f = 0;

    for (const auto& e : arr.edge_handles()) {
        std::list<Point> arc_points;


        arc_points.push_back(e->source()->point());
        arc_points.push_back(e->target()->point());


        int vertex_count = arc_points.size();

        //collect vertices
        double* x = new double[vertex_count];
        double* y = new double[vertex_count];

        int i = 0;
        for (Point p : arc_points) {
            x[i] = to_double(p.x());
            y[i] = to_double(p.y());
            i++;
        }

        //create polygon object
        SHPObject* shape = SHPCreateSimpleObject(SHPT_ARC, vertex_count, x, y, nullptr);

        //write shape into file
        int shape_id = SHPWriteObject(shapefile, -1, shape);

        //set field to face ID
        DBFWriteIntegerAttribute(dbfile, shape_id, field_face_id, f++);

        //memory management
        SHPDestroyObject(shape);
        delete[] x;
        delete[] y;

    }

    DBFClose(dbfile);
    SHPClose(shapefile);



    //create a separate Shapefile for the faces
    //create handle
    shapefile = SHPCreate((path + "_faces").c_str(), SHPT_POLYGON);
    dbfile = DBFCreate((path + "_faces").c_str());

    //create field for arc ID
    field_face_id = DBFAddField(dbfile, "ID", FTInteger, 5, 0);

    f = 0;

    for (const auto& face : arr.face_handles()) {
        std::list<Point> arc_points;

        //ignore outer face
        if (!face->has_outer_ccb()) {
            continue;
        }

        std::vector<Arrangement::Halfedge_const_handle> inner_halfedges;

        Arrangement::Ccb_halfedge_const_circulator first = face->outer_ccb();
        auto curr = first;
        do {
            Arrangement::Halfedge_const_handle e = curr;
            inner_halfedges.push_back(e);
            curr++;
        } while (curr != first);

        for (auto e : inner_halfedges) {

            arc_points.push_back(e->source()->point());
            arc_points.push_back(e->target()->point());

        }

        int vertex_count = arc_points.size();

        //collect vertices
        double* x = new double[vertex_count];
        double* y = new double[vertex_count];

        int i = 0;
        for (Point p : arc_points) {
            x[i] = to_double(p.x());
            y[i] = to_double(p.y());
            i++;
        }

        //create polygon object
        SHPObject* shape = SHPCreateSimpleObject(SHPT_POLYGON, vertex_count, x, y, nullptr);

        //write shape into file
        int shape_id = SHPWriteObject(shapefile, -1, shape);

        //set field to face ID
        DBFWriteIntegerAttribute(dbfile, shape_id, field_face_id, face->data());

        //memory management
        SHPDestroyObject(shape);
        delete[] x;
        delete[] y;

    }

    DBFClose(dbfile);
    SHPClose(shapefile);
}