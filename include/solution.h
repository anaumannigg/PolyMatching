#ifndef _solution_included_
#define _solution_included_

#include "cgal_includes.h"

//class to store a many to many matching as match indices per polygon of both sets
//unmatched polygons get a negative matching index
class Solution {
	std::vector < std::pair<std::vector<int>, std::vector<int>>> matching;
	//assigns every osm polygon an index of the group it belongs to
	std::vector<int> osm_match_index;
	//assigns every osm polygon a match weight corresponding to the weight of the selected edge in E3 by the ILP
	std::vector<double> osm_match_weight;
	//assigns every atkis polygon an index of the group it belongs to
	std::vector<int> atkis_match_index;
	//assigns every atkis polygon a match weight corresponding to the weight of the selected edge in E3 by the ILP
	std::vector<double> atkis_match_weight;
	//remember amount of matches
	int match_count;
	//remember current match ids
	std::pair<int, int> match_ids;
	//remember target value of solution
	double target_value;

public:
    //standard constructor, initializes empty solution
    Solution();

	//init solution instance
	Solution(int num_osm_polys,int num_atkis_polys);

	//init solution instance with specific beginning IDs for matchings, to work on subsets and merge solution afterwards
	Solution(int num_osm_polys, int num_atkis_polys, std::pair<int, int> match_id_limits);

	void addMatch(std::vector<int> set_osm, std::vector<int> set_atkis, double match_weight);

	void insert(Solution sol);
	void insert(Solution sol, std::vector<int> osm_lookup, std::vector<int> atkis_lookup);
	//completes the matching in the way that every non-matched vertex, a negative increasing number is assigned, s.t. they can be distinguished from successfull matches
	void completeMatching();

	int getMatchCount();

	std::pair<int, int> getMatchIDs();

	std::vector<int> getMatchIndices(bool map);

	std::vector<double> getMatchWeights(bool map);

	std::vector < std::pair<std::vector<int>, std::vector<int>>> getMatching();

	double getTargetValue();
};

#endif