#include "../include/solution.h"

Solution::Solution() {
    this->matching = std::vector < std::pair<std::vector<int>, std::vector<int>>>();

    //init vectors to remember which match each polygon belongs to (-1 means unmatched)
    this->osm_match_index = std::vector<int>();
    this->osm_match_weight = std::vector<double>();
    this->atkis_match_index = std::vector<int>();
    this->atkis_match_weight = std::vector<double>();

    this->match_count = 0;

    this->match_ids = { -1,0 };

    this->target_value = 0.0;

}

Solution::Solution(int num_osm_polys,int num_atkis_polys) {
	this->matching = std::vector < std::pair<std::vector<int>, std::vector<int>>>();

	//init vectors to remember which match each polygon belongs to (-1 means unmatched)
	this->osm_match_index = std::vector<int>(); this->osm_match_index.resize(num_osm_polys, -1);
	this->osm_match_weight = std::vector<double>(); this->osm_match_weight.resize(num_osm_polys, 0.0);
	this->atkis_match_index = std::vector<int>(); this->atkis_match_index.resize(num_atkis_polys, -1);
	this->atkis_match_weight = std::vector<double>(); this->atkis_match_weight.resize(num_atkis_polys, 0.0);

	this->match_count = 0;

	this->match_ids = { -1,0 };

	this->target_value = 0.0;
}

Solution::Solution(int num_osm_polys, int num_atkis_polys, std::pair<int,int> match_id_limits) {
	this->matching = std::vector < std::pair<std::vector<int>, std::vector<int>>>();

	//init vectors to remember which match each polygon belongs to (-1 means unmatched)
	this->osm_match_index = std::vector<int>(); this->osm_match_index.resize(num_osm_polys, -1);
	this->osm_match_weight = std::vector<double>(); this->osm_match_weight.resize(num_osm_polys, 0.0);
	this->atkis_match_index = std::vector<int>(); this->atkis_match_index.resize(num_atkis_polys, -1);
	this->atkis_match_weight = std::vector<double>(); this->atkis_match_weight.resize(num_atkis_polys, 0.0);

	this->match_count = 0;

	this->match_ids = match_id_limits;

	this->target_value = 0.0;
}

void Solution::addMatch(std::vector<int> set_osm, std::vector<int> set_atkis, double match_weight) {
	//only add the match, if both sets are non-empty, unmatched polygons should get negative match indices via "completeMatching" in the end
	if (set_osm.size() == 0 || set_atkis.size() == 0) return;

	std::pair<std::vector<int>, std::vector<int>> pair(set_osm, set_atkis);
	this->matching.push_back(pair);
	//remember match index for each poly
	for (const auto& p : set_osm) {
		this->osm_match_index[p] = this->match_ids.second;
		this->osm_match_weight[p] = match_weight;
	}
	for (const auto& p : set_atkis) {
		this->atkis_match_index[p] = this->match_ids.second;
		this->atkis_match_weight[p] = match_weight;
	}
	this->match_ids.second++;
	this->match_count++;

	this->target_value += match_weight;
}

void Solution::insert(Solution sol) {
	auto sol_matching = sol.getMatching();

	//for each match in the argument, add the match to the solution instance
	int i = 0;
	for (const auto& m : sol_matching) {
		std::vector<int> set_osm = m.first;
		std::vector<int> set_atkis = m.second;
		this->addMatch(set_osm, set_atkis, sol.getMatchWeights(0)[set_osm[0]]);
	}

}

void Solution::insert(Solution sol, std::vector<int> osm_lookup, std::vector<int> atkis_lookup) {
	auto sol_matching = sol.getMatching();

	//for each match in the argument, add the match to the solution instance
	int i = 0;
	for (const auto& m : sol_matching) {
		std::vector<int> set_osm, set_atkis;
		for (const auto& i : m.first) set_osm.push_back(osm_lookup[i]);
		for (const auto& i : m.second) set_atkis.push_back(atkis_lookup[i]);
		this->addMatch(set_osm, set_atkis, sol.getMatchWeights(0)[m.first[0]]);
	}

}

void Solution::completeMatching() {
	int match_id = this->match_ids.first;
	for (int i = 0; i < this->osm_match_index.size(); i++) {
		if (this->osm_match_index[i] == -1) {
			//set match ID for unmatched polygon to negative match ID, weight does not need to be updated, as it is 0.0 as default
			this->osm_match_index[i] = match_id--;
			this->match_count++;
		}
	}
	for (int i = 0; i < this->atkis_match_index.size(); i++) {
		if (this->atkis_match_index[i] == -1) {
			//set match ID for unmatched polygon to negative match ID, weight does not need to be updated, as it is 0.0 as default
			this->atkis_match_index[i] = match_id--;
			this->match_count++;
		}
	}

	this->match_ids.first = match_id;
}

int Solution::getMatchCount() {
	return this->match_count;
}

std::pair<int, int> Solution::getMatchIDs() {
	return this->match_ids;
}

//get match indices of the map indicated by the argument
std::vector<int> Solution::getMatchIndices(bool map){
	if (!map) return this->osm_match_index;
	else return this->atkis_match_index;
}

std::vector<double> Solution::getMatchWeights(bool map) {
	if (!map) return this->osm_match_weight;
	else return this->atkis_match_weight;
}

std::vector < std::pair<std::vector<int>, std::vector<int>>> Solution::getMatching() {
	return this->matching;
}

double Solution::getTargetValue() {
	return this->target_value;
}


