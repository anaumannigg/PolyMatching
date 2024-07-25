#include "../include/localization.h"

//LOCALIZATION CLASS

Localization::Localization(std::vector<Polygon_2> polys) {
    //init r-tree
    this->rtree = bgi::rtree< Value_rtree, bgi::quadratic<16>>();

    //init polygon and alpha threshold storage
    this->polys = polys;

    for (int i = 0; i < polys.size(); i++) {

        //input the bounding boxes of the polygons, connected with their index, into the r-tree
        auto bbox = polys[i].bbox();
        //get lowerleft and upperright point of the bounding box
        Point_rtree ll = Point_rtree(bbox.xmin(), bbox.ymin());
        Point_rtree ur = Point_rtree(bbox.xmax(), bbox.ymax());
        //create box for rtree
        Box_rtree b(ll, ur);
        //insert polygon and original index into r-tree
        rtree.insert(std::make_pair(b, i));


    }
}

//locates the input polygon p within the r-tree of the instance and writes neighbors into given vector
//returns true if p has intersecting neighbors within the r tree, else false
bool Localization::get_neighbors(Polygon_2 p, std::vector<Polygon_2>* neighbors) const {
    //compute bounding box of circular arc
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);

    for (const auto& b : query_result) {
        //remember every polygon 
        Polygon_2 query_poly = this->polys[b.second];
        if (CGAL::do_intersect(query_poly, p)) {
            neighbors->push_back(query_poly);
        }
    }

    if (neighbors->size() > 0) return true;
    else return false;
}

//locates the input polygon p within the r-tree of the instance and writes neighbor indices into given vector
//returns true if p has intersecting neighbors within the r tree, else false
bool Localization::get_neighbors(Polygon_2 p, std::vector<int>* neighbors) const {
    //compute bounding box of circular arc
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);

    for (const auto& b : query_result) {
        //remember every polygon 
        Polygon_2 query_poly = this->polys[b.second];
        if (CGAL::do_intersect(query_poly, p)) {
            neighbors->push_back(b.second);
        }
    }
    //sort output vector
    std::sort(neighbors->begin(), neighbors->end());
    //make output vector unique
    neighbors->erase(std::unique(neighbors->begin(), neighbors->end()), neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

//locates the input point and writes intersected polygons into neighbors
bool Localization::get_neighbors(Point p, std::vector<int>* neighbors) const {
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);

    for (const auto& b : query_result) {
        //remember every polygon 
        Polygon_2 query_poly = this->polys[b.second];
        if (!query_poly.has_on_unbounded_side(p)) {
            neighbors->push_back(b.second);
        }
    }
    //sort output vector
    std::sort(neighbors->begin(), neighbors->end());
    //make output vector unique
    neighbors->erase(std::unique(neighbors->begin(), neighbors->end()), neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

//locates the input segment and writes intersected polygons into neighbors
//returns true if p has intersecting neighbors, else false
bool Localization::get_neighbors(Segment p, std::vector<int>* neighbors) const {
    //if the segment represents a single point, return
    if (p.source() == p.target()) return this->get_neighbors(p.source(),neighbors);

    //build poly from segment
    Polygon_2 p_poly;
    p_poly.push_back(p.source()); p_poly.push_back(p.target());

    //compute bounding box of circular arc
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);

    for (const auto& b : query_result) {
        //retreive query polygon
        Polygon_2 query_poly = this->polys[b.second];

        //make sure the input segment is not one segment on the boundary of the polygon, as then the intersection check will crash
        bool is_boundary_edge = false;
        for (const auto& e : query_poly.edges()) {
            if (e.source() == p.source() && e.target() == p.target() || e.target() == p.source() && e.source() == p.target()) {
                is_boundary_edge = true;
                break;
            }
        }
        if (is_boundary_edge) {
            //collision found
            neighbors->push_back(b.second);
        }
        //input segment is not an edge of the polygon, check for intersection
        else {
            Polygon_set_2 query_poly_set;
            query_poly_set.insert(query_poly);
            if (!query_poly_set.do_intersect(p_poly)) {
                neighbors->push_back(b.second);
            }
        }
    }
    //sort output vector
    std::sort(neighbors->begin(), neighbors->end());
    //make output vector unique
    neighbors->erase(std::unique(neighbors->begin(), neighbors->end()), neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

//neighbors function with constraint: an intersecting polygon is only returned as a neighbor, iff union / smaller_poly_area is >= the argument threshold
bool Localization::get_neighbors(Polygon_2 p, std::vector<int>* neighbors, double neighboring_threshold) const {
    //compute bounding box of circular arc
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);

    for (const auto& b : query_result) {
        //remember every polygon 
        Polygon_2 query_poly = this->polys[b.second];

        
        if (CGAL::do_intersect(query_poly, p)) {
            //polys do intersect,but further requirements have to be met
            //compute intersection area
            std::list<Polygon_wh_2> inter;
            CGAL::intersection(query_poly, p, std::back_inserter(inter));
            double inter_area = 0.0;
            for (const auto& inter_p : inter) {
                inter_area += to_double(inter_p.outer_boundary().area());
                for (const auto& hole : inter_p.holes()) {
                    inter_area -= to_double(hole.area());
                }
            }

            //get smaller polygon area
            double smaller_poly_area = std::min(to_double(query_poly.area()), to_double(p.area()));

            if(inter_area/smaller_poly_area > neighboring_threshold) neighbors->push_back(b.second);
        }
    }
    //sort output vector
    std::sort(neighbors->begin(),neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

//neighbors function with constraint: an intersecting polygon is only returned as a neighbor, iff union / smaller_poly_area is >= the argument threshold, speedup via arrangement
bool Localization::get_neighbors(Polygon_2 p, int p_index, bool p_map, std::vector<int>* neighbors, double neighboring_threshold, MapOverlay mo) const {
    //compute bounding box of circular arc
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);


    

    for (const auto& b : query_result) {
        //remember every polygon 
        Polygon_2 query_poly = this->polys[b.second];

        if (CGAL::do_intersect(p, query_poly)) {

            std::vector<int> osm_indices = !p_map ? std::vector<int> {p_index} : std::vector<int> {(int)b.second};
            std::vector<int> atkis_indices = !p_map ? std::vector<int> {(int)b.second} : std::vector<int>{ p_index };

            double inter_area = mo.getIoU(osm_indices, atkis_indices);

            if (inter_area > 0) {

                //get smaller polygon area
                double smaller_poly_area = std::min(to_double(query_poly.area()), to_double(p.area()));

                if (inter_area / smaller_poly_area > neighboring_threshold) neighbors->push_back(b.second);
            }
        }
    }
    //sort output vector
    std::sort(neighbors->begin(), neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

bool Localization::get_neighbors_fully_included(Polygon_2 p, std::vector<int>* neighbors) const {
    //compute bounding box of circular arc
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);

    for (const auto& b : query_result) {
        //remember every polygon, but only if it has not ben considered yet and the query poly lies entirely within it
        Polygon_2 query_poly = this->polys[b.second];

        if (std::find(neighbors->begin(),neighbors->end(),b.second) == neighbors->end()) {
            //polygon has not been considered yet as neighbor
            //check if every point of p lies inside the polygon (more stable and efficient than intersection area check)
            bool all_points_lie_inside = true;
            for (const auto& point : p.vertices()) {
                if (query_poly.has_on_unbounded_side(point)) {
                    all_points_lie_inside = false;
                    break;
                }
            }

            if (all_points_lie_inside) neighbors->push_back(b.second);

        }
    }
    //sort output vector
    std::sort(neighbors->begin(), neighbors->end());
    //make output vector unique
    neighbors->erase(std::unique(neighbors->begin(), neighbors->end()), neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

bool Localization::get_neighbors_majorly_included(Polygon_2 p, std::vector<int>* neighbors, double inclusion_threshold) const {
    //compute bounding box of circular arc
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);

    for (const auto& b : query_result) {
        //remember every polygon, but only if it has not ben considered yet and the query poly lies entirely within it
        Polygon_2 query_poly = this->polys[b.second];

        if (std::find(neighbors->begin(), neighbors->end(), b.second) == neighbors->end()) {
            //polygon has not been considered yet as neighbor
            
            //compute o_tilde for p and b
            std::list<Polygon_wh_2> inter;
            CGAL::intersection(query_poly, p, std::back_inserter(inter));
            double inter_area = 0.0;
            for (const auto& inter_p : inter) {
                inter_area += to_double(inter_p.outer_boundary().area());
                for (const auto& hole : inter_p.holes()) {
                    inter_area -= to_double(hole.area());
                }
            }

            //get smaller polygon area
            double smaller_poly_area = std::min(to_double(query_poly.area()), to_double(p.area()));

            if (inter_area / smaller_poly_area > inclusion_threshold) neighbors->push_back(b.second);
     
        }
    }
    //sort output vector
    std::sort(neighbors->begin(), neighbors->end());
    //make output vector unique
    neighbors->erase(std::unique(neighbors->begin(), neighbors->end()), neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

bool Localization::get_neighbors_with_minimum_intersection_proportion(Polygon_2 p, std::vector<int>* neighbors, double inclusion_threshold) const {
    //compute bounding box of circular arc
    auto poly_bbox = p.bbox();
    Point_rtree ll = Point_rtree(poly_bbox.xmin(), poly_bbox.ymin());
    Point_rtree ur = Point_rtree(poly_bbox.xmax(), poly_bbox.ymax());
    Box_rtree query_box = Box_rtree(ll, ur);

    //query r-tree
    std::vector<Value_rtree> query_result = this->query(query_box);

    for (const auto& b : query_result) {
        //remember every polygon, but only if it has not ben considered yet and the query poly lies entirely within it
        Polygon_2 query_poly = this->polys[b.second];

        if (std::find(neighbors->begin(), neighbors->end(), b.second) == neighbors->end()) {
            //polygon has not been considered yet as neighbor

            //compute o_tilde for p and b
            std::list<Polygon_wh_2> inter;
            CGAL::intersection(query_poly, p, std::back_inserter(inter));
            double inter_area = 0.0;
            for (const auto& inter_p : inter) {
                inter_area += to_double(inter_p.outer_boundary().area());
                for (const auto& hole : inter_p.holes()) {
                    inter_area -= to_double(hole.area());
                }
            }
            //check if at last the necessary threshold of the query poly area is in the intersection
            if (inter_area / to_double(query_poly.area()) > inclusion_threshold) neighbors->push_back(b.second);

        }
    }
    //sort output vector
    std::sort(neighbors->begin(), neighbors->end());
    //make output vector unique
    neighbors->erase(std::unique(neighbors->begin(), neighbors->end()), neighbors->end());
    if (neighbors->size() > 0) return true;
    else return false;
}

std::vector<Value_rtree> Localization::query(Box_rtree query_box) const {
    std::vector<Value_rtree> query_result;
    this->rtree.query(bgi::intersects(query_box),std::back_inserter(query_result));
    return query_result;
}

bool Localization::are_adjacent(Polygon_2 polyA, int polyB_index) const {
    //get neighbors of poly A
    std::vector<int> neighbors;
    //note: here, different get_neighbors modes could be applied in order to adapt overlap checks
    this->get_neighbors(polyA, &neighbors);

    //check, if polyB is in the neighbors list of A
    if (std::find(neighbors.begin(), neighbors.end(), polyB_index) != neighbors.end()) return true;
    else return false;
}

//END LOCALIZATION CLASS

// MAP OVERLAY CLASS

MapOverlay::MapOverlay(std::vector<Polygon_2> osm_polys, std::vector<Polygon_2> atkis_polys) {
    this->arr = Arrangement();

    //insert all segments of polygons into the arrangement
    std::vector<ArrSegment> segments;
    for (const auto& polys : { osm_polys,atkis_polys }) {
        for (const auto& poly : polys) {
            for (const auto& e : poly.edges()) {
                segments.push_back(e);
            }
        }
    }
    CGAL::insert(arr, segments.begin(), segments.end());

    this->areas = std::vector<double>();
    this->areas.resize(arr.number_of_faces(), 0.0);
    this->osm_poly_covered_faces = std::vector<std::vector<int>>(osm_polys.size());
    this->atkis_poly_covered_faces = std::vector<std::vector<int>>(atkis_polys.size());


}

void MapOverlay::assignFaces(std::vector<Polygon_2> osm_polys, std::vector<Polygon_2> atkis_polys, Localization osm_rtree, Localization atkis_rtree) {
    int face_id = 0;

    for (auto fit = this->arr.faces_begin(); fit != this->arr.faces_end(); fit++) {


        //if the face is the outer face or has no outer ccb, do not compute the face area
        if (fit->is_unbounded() || !fit->has_outer_ccb()) {
            face_id++;
            continue;
        }

        //compute face area
        double face_area = 0.0;

        //first compute and add outer ccb
        Polygon_2 outer_ccb;
        Arrangement::Ccb_halfedge_const_circulator e = fit->outer_ccb();
        auto loop = e;
        do {
            auto x = e->source();
            outer_ccb.push_back(e->source()->point());
        } while (++e != loop);
        //make sure that the area can actually be computed, there might be artifacts leading to a self-intersecting outer ccb
        if (outer_ccb.is_simple()) {
            face_area += to_double(outer_ccb.area());

            //write area into vector
            this->areas[face_id] = face_area;


            //check, which polygons overlap the current face
            std::vector<int> osm_neighbors, atkis_neighbors;
            osm_rtree.get_neighbors(outer_ccb, &osm_neighbors);
            atkis_rtree.get_neighbors(outer_ccb, &atkis_neighbors);

            //remember for those polygons, that the facet is included within them 
            for (const auto& osm_n : osm_neighbors) this->osm_poly_covered_faces[osm_n].push_back(face_id);
            for (const auto& atkis_n : atkis_neighbors) this->atkis_poly_covered_faces[atkis_n].push_back(face_id);


            face_id++;
        }
    }

}

std::vector<int> MapOverlay::getCoveredFaces(bool map, int poly_index) const {
    std::vector<std::vector<int>> polys = !map ? this->osm_poly_covered_faces : this->atkis_poly_covered_faces;
    return polys[poly_index];
}

std::vector<int> MapOverlay::getCoveredFaces(bool map, std::vector<int> poly_indices) const {
    std::vector<int> covered_faces;
    for (const auto& poly : poly_indices) {
        std::vector<int> cov_faces_of_poly = getCoveredFaces(map, poly);
        if (cov_faces_of_poly.size() > 0) covered_faces.insert(covered_faces.end(), cov_faces_of_poly.begin(), cov_faces_of_poly.end());
    }
    return covered_faces;
}

double MapOverlay::getIoU(std::vector<int> osm_indices, std::vector<int> atkis_indices) const {
    //get common intersected faces of polygons of both sets
    std::vector<int> osm_faces, atkis_faces;
    for (const auto& osm_i : osm_indices) {
        std::vector<int> covered_faces = this->getCoveredFaces(0, osm_i);
        osm_faces.insert(osm_faces.end(), covered_faces.begin(), covered_faces.end());
    }
    for (const auto& atkis_i : atkis_indices) {
        std::vector<int> covered_faces = this->getCoveredFaces(1, atkis_i);
        atkis_faces.insert(atkis_faces.end(), covered_faces.begin(), covered_faces.end());
    }

    return this->getIoUviaCoveredFaces(osm_faces, atkis_faces);
}

double MapOverlay::getIntersectionArea(std::vector<int> osm_indices, std::vector<int> atkis_indices) const {
    //get common intersected faces of polygons of both sets
    std::vector<int> osm_faces, atkis_faces;
    for (const auto& osm_i : osm_indices) {
        std::vector<int> covered_faces = this->getCoveredFaces(0, osm_i);
        osm_faces.insert(osm_faces.end(), covered_faces.begin(), covered_faces.end());
    }
    for (const auto& atkis_i : atkis_indices) {
        std::vector<int> covered_faces = this->getCoveredFaces(1, atkis_i);
        atkis_faces.insert(atkis_faces.end(), covered_faces.begin(), covered_faces.end());
    }

    return this->getIntersectionviaCoveredFaces(osm_faces, atkis_faces);
}

double MapOverlay::getIntersectionviaCoveredFaces(std::vector<int> osm_faces, std::vector<int> atkis_faces) const {
    //sort both lists and make unique
    std::sort(osm_faces.begin(), osm_faces.end());
    osm_faces.erase(unique(osm_faces.begin(), osm_faces.end()), osm_faces.end());
    std::sort(atkis_faces.begin(), atkis_faces.end());
    atkis_faces.erase(unique(atkis_faces.begin(), atkis_faces.end()), atkis_faces.end());

    //get intersection and union of face indices of both sets
    std::vector<int> inter_faces;
    std::set_intersection(osm_faces.begin(), osm_faces.end(), atkis_faces.begin(), atkis_faces.end(), std::back_inserter(inter_faces));


    double inter_area = 0.0;
    for (const auto& inter_face : inter_faces) inter_area += this->areas[inter_face];
    return inter_area;
}

double MapOverlay::getUnionviaCoveredFaces(std::vector<int> osm_faces, std::vector<int> atkis_faces) const {
    //sort both lists and make unique
    std::sort(osm_faces.begin(), osm_faces.end());
    osm_faces.erase(unique(osm_faces.begin(), osm_faces.end()), osm_faces.end());
    std::sort(atkis_faces.begin(), atkis_faces.end());
    atkis_faces.erase(unique(atkis_faces.begin(), atkis_faces.end()), atkis_faces.end());

    //get intersection and union of face indices of both sets
    std::vector<int> union_faces;
    union_faces.insert(union_faces.end(), osm_faces.begin(), osm_faces.end());
    union_faces.insert(union_faces.end(), atkis_faces.begin(), atkis_faces.end());
    sort(union_faces.begin(), union_faces.end());
    union_faces.erase(unique(union_faces.begin(), union_faces.end()), union_faces.end());


    double union_area = 0.0;
    for (const auto& union_face : union_faces) union_area += this->areas[union_face];
    return union_area;
}

double MapOverlay::getIoUviaCoveredFaces(std::vector<int> osm_faces, std::vector<int> atkis_faces) const {
    //sort both lists and make unique
    std::sort(osm_faces.begin(), osm_faces.end());
    osm_faces.erase(unique(osm_faces.begin(), osm_faces.end()), osm_faces.end());
    std::sort(atkis_faces.begin(), atkis_faces.end());
    atkis_faces.erase(unique(atkis_faces.begin(), atkis_faces.end()), atkis_faces.end());

    //get intersection and union of face indices of both sets
    std::vector<int> inter_faces, union_faces;
    std::set_intersection(osm_faces.begin(), osm_faces.end(), atkis_faces.begin(), atkis_faces.end(), std::back_inserter(inter_faces));
    union_faces.insert(union_faces.end(), osm_faces.begin(), osm_faces.end());
    union_faces.insert(union_faces.end(), atkis_faces.begin(), atkis_faces.end());
    sort(union_faces.begin(), union_faces.end());
    union_faces.erase(unique(union_faces.begin(), union_faces.end()), union_faces.end());


    double inter_area = 0.0, union_area = 0.0;
    for (const auto& inter_face : inter_faces) inter_area += this->areas[inter_face];

    for (const auto& union_face : union_faces) union_area += this->areas[union_face];
    return inter_area / union_area;
}

double MapOverlay::getArea(bool map, std::vector<int> indices) const {
    std::vector<int> faces;
    for (const auto& i : indices) {
        std::vector<int> covered_faces = this->getCoveredFaces(map, i);
        faces.insert(faces.end(), covered_faces.begin(), covered_faces.end());
    }

    return this->getAreaofCoveredFaces(faces);

}

double MapOverlay::getAreaofCoveredFaces(std::vector<int> covered_faces) const {
    double area = 0.0;
    for (const auto& f : covered_faces) area += this->areas[f];
    return area;
}

bool MapOverlay::doOverlap(std::vector<int> osm_indices, std::vector<int> atkis_indices, double epsilon) const {
    //get common intersected faces of polygons of both sets
    std::vector<int> osm_faces, atkis_faces;

    for (const auto& osm_i : osm_indices) {
        std::vector<int> covered_faces = this->getCoveredFaces(0, osm_i);
        osm_faces.insert(osm_faces.end(), covered_faces.begin(), covered_faces.end());
    }
    for (const auto& atkis_i : atkis_indices) {
        std::vector<int> covered_faces = this->getCoveredFaces(1, atkis_i);
        atkis_faces.insert(atkis_faces.end(), covered_faces.begin(), covered_faces.end());
    }


    return this->doOverlapviaCoveredFaces(osm_faces, atkis_faces, epsilon);

}

bool MapOverlay::doOverlapviaCoveredFaces(std::vector<int> osm_covered_faces, std::vector<int> atkis_covered_faces, double epsilon) const {
    double inter_area = this->getIntersectionviaCoveredFaces(osm_covered_faces, atkis_covered_faces);

    double min_area = std::min(this->getAreaofCoveredFaces(osm_covered_faces), this->getAreaofCoveredFaces(atkis_covered_faces));

    return inter_area / min_area > epsilon;
}

// END MAP OVERLAY CLASS

//helper struct for the comparison of segments to keep a unique set
struct SegmentComparator {
    bool operator()(const Segment& s1, const Segment& s2) const {
        // Compare the x values of the source points first
        if (s1.source().x() != s2.source().x()) {
            return s1.source().x() < s2.source().x();
        }
        // If x values are equal, compare the y values of the source points
        if (s1.source().y() != s2.source().y()) {
            return s1.source().y() < s2.source().y();
        }
        // If the source points are the same, compare the target points
        if (s1.target().x() != s2.target().x()) {
            return s1.target().x() < s2.target().x();
        }
        return s1.target().y() < s2.target().y();
    }
};

// Function to compute the union of two vectors of polygons
std::vector<Polygon_2> merge(const std::vector<Polygon_2>& polys1, const std::vector<Polygon_2>& polys2) {
    //create a set of edges of both sets of polygons
    std::set<Segment, SegmentComparator> segments;
    for (const auto &polys: {polys1, polys2}) {
    for (const auto &p: polys) {
            for(const auto& e : p.edges()) {
                segments.insert(e);
            }
        }
    }

    //create an arrangement using all segments
    Arrangement arr;
    CGAL::insert(arr,segments.begin(),segments.end());

    //label the faces of the arrangement with their ids
    int f_id = 0;
    for(auto f = arr.faces_begin(); f!= arr.faces_end(); f++) {
        f->set_data(f_id++);
    }

    //create an r-tree with all polygons to check, if a face lies on a polygon
    std::vector<Polygon_2> all_polys;
    for(auto& polys : {polys1,polys2}) all_polys.insert(all_polys.end(),polys.begin(),polys.end());
    Localization idx(all_polys);


    //collect all faces, that lie within at least one polygon of the two sets
    std::vector<Arrangement::Face_const_handle> poly_faces;
    std::vector<int> poly_face_ids;
    for(auto f = arr.faces_begin(); f!= arr.faces_end(); f++) {
        //get centroid of the face
        if(f->has_outer_ccb()) {
            Polygon_2 f_p;

            auto e = f->outer_ccb();
            auto e_loop = e;
            do {
                f_p.push_back(e->source()->point());
            }while (++e != e_loop);


            //create a point sample inside the polygon for r-tree query
            //(from a vertex, turn inside the polygon and sample a point)
            Point sample;
            auto v_c = f_p.vertices_circulator();
            do {
                Point source = *v_c;
                Point target_left = *(v_c-1);
                Point target_right = *(v_c+1);
                Vector dir_left(source,target_left), dir_right(source,target_right);
                double angle = get_angle(target_right,target_left);
                Vector inside = rotate(dir_right,angle/2);
                double eps = 1e-10;
                sample = source + eps * normalize(inside);

                //cout << "POINT((" << sample << "))" << endl;

                v_c++;

            } while(!f_p.has_on_bounded_side(sample));

            //if the sampled point lies at least within one polygon, remember the face
            std::vector<int> neighbors;
            if(idx.get_neighbors(sample,&neighbors)) {
                poly_faces.push_back(f);
                poly_face_ids.push_back(f->data());
            }


        }
    }

    //sort poly face id vector for quick check if a face is a poly face
    std::sort(poly_face_ids.begin(),poly_face_ids.end());

    //initialize a vector to keep track of visited faces
    std::vector<bool> visited_faces;
    visited_faces.resize(arr.number_of_faces(),false);

    //find hole faces, that are completely surrounded by poly faces
    //as we only want to find the outer boundary of the union, holes should be neglected
    for(auto f = arr.faces_begin(); f != arr.faces_end(); f++) {
        //check if it is a bounded face and if it has not been labeled a polygon face already
        if(f->has_outer_ccb() && !std::binary_search(poly_face_ids.begin(),poly_face_ids.end(),f->data())) {
            bool is_hole = true;

            auto e = f->outer_ccb();
            auto e_loop = e;
            do {
                if(!std::binary_search(poly_face_ids.begin(),poly_face_ids.end(),e->twin()->face()->data())) {
                    is_hole = false;
                    break;
                }
            }while(++e!=e_loop);

            //mark hole as poly face, such that it does not form a boundary in the trace procedure
            if(is_hole) {
                poly_faces.push_back(f);
                poly_face_ids.push_back(f->data());
            }

        }
    }
    std::sort(poly_face_ids.begin(),poly_face_ids.end());

    //initialize an empty vector for the merged polygons to return later on
    std::vector<Polygon_2> merged_polys;


    //start out with a poly face and trace the boundary as long as finding other selected faces and not arriving at the starting edge again
    for (const auto& p_face : poly_faces) {
        if (!visited_faces[p_face->data()] && p_face->number_of_outer_ccbs() == 1) {

            //mark face as selected
            visited_faces[p_face->data()] = true;

            //check if a selected face, which lies on the boundary of a connected subset is found
            bool is_on_boundary = false;
            Arrangement::Halfedge_const_handle curr;
            auto e_it = p_face->outer_ccb();
            auto e_loop = e_it;
            do {
                Arrangement::Halfedge_const_handle e = e_it;
                if (!std::binary_search(poly_face_ids.begin(),poly_face_ids.end(), e->twin()->face()->data())) {
                    curr = e;
                    is_on_boundary = true;
                    break;
                }
            } while (++e_it != e_loop);


            //if edge on boundary is found, begin trace
            if (is_on_boundary) {

                //prepare polygon to push back into solution
                Polygon_2 poly;

                auto loop = curr;
                do {
                    //push back next point
                    poly.push_back(curr->source()->point());
                    //continue trace of the boundary
                    curr = curr->next();
                    while (std::binary_search(poly_face_ids.begin(), poly_face_ids.end(), curr->twin()->face()->data())) {
                        curr = curr->twin()->next();
                        //mark face as visited
                        visited_faces[curr->face()->data()] = true;
                    }

                    //cout << "LINESTRING((" << std::fixed << to_double(curr->source()->point().x()) << " " << to_double(curr->source()->point().y()) << " , " << to_double(curr->target()->point().x()) << " " << to_double(curr->target()->point().y()) << "))" << endl;

                } while (curr != loop);

                //add polygon to list
                merged_polys.push_back(poly);

            }
        }


    }

    return merged_polys;
}