#ifndef _localization_included_
#define _localization_included_

#include "cgal_includes.h"
#include "shapefile_io_operations.h"

//for r-tree
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<float, 2, bg::cs::cartesian> Point_rtree;
typedef bg::model::box<Point_rtree> Box_rtree;
typedef std::pair<Box_rtree, unsigned> Value_rtree;

class Localization;
class MapOverlay;

//a class storing an r-tree for efficient spacial queries
class Localization {
	//r tree for fast collection of polygons to check
	bgi::rtree<Value_rtree, bgi::quadratic<16>> rtree;
	//polys
	std::vector<Polygon_2> polys;

public:
	Localization(std::vector<Polygon_2> polys);
	bool get_neighbors(Polygon_2 p, std::vector<Polygon_2>* neighbors) const;
    bool get_neighbors(Polygon_2 p, std::vector<int>* neighbors) const;
    bool get_neighbors(Point p, std::vector<int>* neighbors) const;
    bool get_neighbors(Segment p, std::vector<int>* neighbors) const;
    bool get_neighbors(Polygon_2 p, std::vector<int>* neighbors, double neighboring_threshold) const;
    bool get_neighbors(Polygon_2 p, int p_index, bool p_map, std::vector<int>* neighbors, double neighboring_threshold, MapOverlay mo) const;
	//writes all neighbors into the neighbors vector, for which the input polygon p lies entirely within them, returns true if neighbors have been found, else false
    bool get_neighbors_fully_included(Polygon_2 p, std::vector<int>* neighbors) const;
	//returns all neighbors n, where I(p,n) / min(area(p),area(n)) > inclusion_threshold
    bool get_neighbors_majorly_included(Polygon_2 p, std::vector<int>* neighbors, double inclusion_threshold) const;
	//returns all neighbors n, where at least inclusion_threshold [0,1] of the area of n is included in query poly p
    bool get_neighbors_with_minimum_intersection_proportion(Polygon_2 p, std::vector<int>* neighbors, double inclusion_threshold) const;
    bool are_adjacent(Polygon_2 polyA, int polyB_index) const;
    std::vector<Value_rtree> query(Box_rtree query_box) const;

private:
	std::vector<Value_rtree> query(Box_rtree query_box);
};

//class to create and store and arrangement of all polygons to reduce computation effort of area computations in intersection over union
class MapOverlay {
	Arrangement arr;
	std::vector<double> areas;
	std::vector<std::vector<int>> osm_poly_covered_faces, atkis_poly_covered_faces;

public:
    //computes an arrangement of all edges in both sets of polgons. the area of every face gets computed
    //such that intersection- and union areas can then be computed via simple lookups instead of geometric operations.
	MapOverlay(std::vector<Polygon_2> osm_polys, std::vector<Polygon_2> atkis_polys);
	void assignFaces(std::vector<Polygon_2> osm_polys, std::vector<Polygon_2> atkis_polys, Localization osm_rtree, Localization atkis_rtree);
	std::vector<int> getCoveredFaces(bool map, int poly_index) const;
    std::vector<int> getCoveredFaces(bool map, std::vector<int> poly_indices) const;
    double getIoU(std::vector<int> osm_indices, std::vector<int> atkis_indices) const;
    double getIntersectionArea(std::vector<int> osm_indices, std::vector<int> atkis_indices) const;
    double getIntersectionviaCoveredFaces(std::vector<int> osm_covered_faces, std::vector<int> atkis_covered_faces) const;
    double getUnionviaCoveredFaces(std::vector<int> osm_covered_faces, std::vector<int> atkis_covered_faces) const;
    double getIoUviaCoveredFaces(std::vector<int> osm_covered_faces, std::vector<int> atkis_covered_faces) const;
    double getArea(bool map, std::vector<int> indices) const;
    double getAreaofCoveredFaces(std::vector<int> covered_faces) const;
    bool doOverlap(std::vector<int> osm_indices, std::vector<int> atkis_indices, double epsilon) const;
    bool doOverlapviaCoveredFaces(std::vector<int> osm_covered_faces, std::vector<int> atkis_covered_faces, double epsilon) const;
};

std::vector<Polygon_2> merge(const std::vector<Polygon_2>& polys1, const std::vector<Polygon_2>& polys2);

#endif