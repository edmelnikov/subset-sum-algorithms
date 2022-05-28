/*
Implementation of Bringmann's algorithm described in "A Near-Linear Pseudopolynomial Time Algorithm for Subset Sum" (2016)
*/
#include <iostream>
#include <stdio.h>
#include <complex.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <time.h>
#include <utility>
#include <chrono>
#include <random>
#include <unordered_map>

#include "fftw3.h"
#include "minkowski_sum.h"
#include "n_choose_k.h"
#include "bellman_ssum.h"

unsigned int seed = std::chrono::steady_clock::now().time_since_epoch().count();
std::mt19937 r_eng(seed);

/* Randomly partitions given "set" into "num_parts" subsets */
std::unordered_map<int, std::vector<int>> rand_part(const std::vector<int>& set, int num_parts) {
	std::unordered_map<int, std::vector<int>> set_partitioned; // hash table
	std::uniform_int_distribution<int> distr(0, num_parts - 1);
	
	int color;
	for (auto element : set) { 
		//color = rand() % num_parts; // assign a random color
		color = distr(r_eng); 
		set_partitioned[color].push_back(element); 
	}
	return set_partitioned;
}

/* 
Algorithm 1. Computes a set of subset sums bounded by "target" of a "set" such that
any subset of a "set" with cardinality less than or equal to "k" and its sum less than or equal to "target"
belongs to a set of all subset sums bounded by "target" of a "set" with probability greater than 1 - "delta"
*/
std::vector<int> ColorCoding(const std::vector<int>& set, int target, int k, double delta) {
	/* Handle trivial cases */
	if (k == 0) return std::vector<int>({});
	if (k == 1 || set.size() == 0 || set.size() == 1) return set;
	if (set.size() == 2) return std::vector<int>({ set.front(), set.back(), set.front() + set.back() });
	if (k >= set.size()) {
		return get_ssums_dyn_table(set, target);
	}

	// int num_iter = ceil((log((1 / delta)) / log((double)4/3)));
	int num_iter = ceil(log2(1 / delta));
	std::vector<int> mink_prod_united, mink_prod_united_old;
	for (int i = 0; i < num_iter; i++) { // iterate over iterations
		std::unordered_map<int, std::vector<int>> set_partitioned = std::move(rand_part(set, k*k)); // randomly partition into k^2 sets
		
		/* Calculate minkowsky sum across all subsets of a partition */
		std::vector<int> mink_prod = {};
		for (auto subset : set_partitioned) { // iterate over map's key-value pairs and apply minkowski sum to all subsets
			//mink_prod = std::move(minkowski_add(mink_prod, subset.second, target));
			//mink_prod = std::move(minkowski_sum_mpir(mink_prod, subset.second, target));
			mink_prod = std::move(minkowski_sum_fft(mink_prod, subset.second, target));
		}
		
		/* Unite */
		mink_prod_united_old = std::move(mink_prod_united); // move elements to the buffer "mink_prod_united_old" in order to record the result to "mink_prod_united"
		std::set_union(
			mink_prod_united_old.begin(),
			mink_prod_united_old.end(),
			mink_prod.begin(),
			mink_prod.end(),
			std::back_inserter(mink_prod_united)
		); 
	}
	return mink_prod_united;
}

/*
Algorithm 3. Given a "layer"-layer (l-layer) "set", computes its set of all subset sums bounded by "target" containing 
any sum with probability at least 1 - "delta"
*/
std::vector<int> ColorCodingLayer(const std::vector<int>& set, int target, int layer, double delta) {
	/* Handle trivial cases */
	if (set.size() == 0 || set.size() == 1) return set;
	if (set.size() == 2) return std::vector<int>({ set.front(), set.back(), set.front() + set.back() });

	/* Calculate log(ℓ/δ)*/
	double log_layer_delta = log2((double)layer / delta);

	/* Check if ℓ < log(ℓ/δ) */
	if (layer < log_layer_delta) {
		return ColorCoding(set, target, layer, delta);
	}
	 
	/* Calculate number of subsets */
	int num_part = pow2_bs(ceil(log2(layer / log_layer_delta)));
	
	/* Calculate "k" parameter of a subset */
	int subset_k = 6 * log_layer_delta;

	/* Randomly partition given "set" into "num_part" subsets */
	std::unordered_map<int, std::vector<int>> set_partitioned = std::move(rand_part(set, num_part));

	/* Calculte subset sums of each subset of the partition with color coding */
	std::vector<std::vector<int>> set_partitioned_ssums(set_partitioned.size(), std::vector<int>()); // vector of subset sums of each subset of the partition
	int set_partitioned_ssums_ind = 0;	
	for (auto subset : set_partitioned) { // iterate over subsets of a partitioned set

		/* Calculate subset target and check if it's larger than the original target */
		int subset_target = ceil(2 * subset_k * target / layer); // target of a subset
		if (subset_target > target) subset_target = target;

		/* Apply ColorCoding algorthm to the subset with the specified parameters */
		set_partitioned_ssums[set_partitioned_ssums_ind] = std::move(ColorCoding(
			subset.second, // use .second to access the subset itself since subset variable is a std::pair<int, int>
			subset_target,
			subset_k,
			delta/layer
			));
		set_partitioned_ssums_ind++;
	}

	/* apply minkowski sum to all "set_partitioned_ssums" in a binary-tree-like way */
	std::vector<std::vector<int>> curr_round = std::move(set_partitioned_ssums); 
	std::vector<std::vector<int>> next_round;
	for (int round = 0; round < log2(num_part); round++) { // iterate over levels of a tree from bottom to top (round = depth - level)
		for (int j = 0; j < curr_round.size(); j += 2) {
			if (j + 1 < curr_round.size()) {
				//next_round.push_back(minkowski_add(curr_round[j], curr_round[j+1], pow(2, round + 1)*2*6*log_layer_delta*target/layer));
				//next_round.push_back(std::move(minkowski_sum_mpir(curr_round[j], curr_round[j + 1], pow2_bs(round + 1) * 2 * 6 * log_layer_delta * target / layer)));
				next_round.push_back(std::move(minkowski_sum_fft(curr_round[j], curr_round[j + 1], pow2_bs(round + 1) * 2 * 6 * log_layer_delta * target / layer)));
			}
			else {
				next_round.push_back(curr_round[j]);
			}
		}
		curr_round = std::move(next_round);
	}

	/* intersect with {0, ... , target}*/
	std::vector<int> res;
	if (curr_round.size() != 0) {
		res.reserve(curr_round[0].size());
		for (auto el : curr_round[0]) {
			if (el <= target) {
				res.push_back(el);
			}
		}
	}
	return res;
}

/* 
Algorithm 2. Computes a set of subset sums of a given "set" containing any sum with probability at least 1 - "delta"
*/
std::pair<bool, double> bringmann_ssum(const std::vector<int>& set, int target, double delta) {
	auto start_time = std::chrono::high_resolution_clock::now(); // time measurement

	/* Calculate number of layers */
	int num_layers = ceil(log2(set.size())); // calculate number of layers
	
	std::vector<int> layer_set; // a vector containing a set of the current layer 
	layer_set.reserve(set.size()); // reserve enough space
	int set_ind = 0; // index of a set's element
	double right_bound = (double)target / pow2_bs(num_layers - 1); // right bound of the current layer
	std::vector<int> layer_ssums; // all subset sums of the current layer
	std::vector<int> all_ssums = {}; // subset sums of all layers 

	/* iterate over layers backwards, adding elements from the set in the process */
	for (int layer_ind = num_layers - 1; layer_ind >= 0; layer_ind--) {
		// note: this code ignores elements higher than target 
		while (set[set_ind] <= right_bound && set_ind < set.size()) { // check if current set's item fits into the layer
			layer_set.push_back(set[set_ind]); // if it does, add it to the current layer
			set_ind++; // and increment the index of a set
		}
		right_bound = (double)target / pow2_bs(layer_ind - 1); // recalculate the right bound once all the items have been added to the current layer
		layer_ssums = std::move(ColorCodingLayer(layer_set, target, pow2_bs(layer_ind + 1), delta / num_layers));
		
		//all_ssums = std::move(minkowski_add(all_ssums, layer_ssums, target)); // note: minkowsky sum always returns a sorted set
		//all_ssums = std::move(minkowski_sum_mpir(all_ssums, layer_ssums, target));
		all_ssums = std::move(minkowski_sum_fft(all_ssums, layer_ssums, target));
		layer_set.clear();
	}

	bool solution = 0;
	if (std::find(all_ssums.begin(), all_ssums.end(), target) != all_ssums.end()) {
		solution = 1;
	}

	auto end_time = std::chrono::high_resolution_clock::now(); // time measurement
	std::chrono::duration<double, std::milli> running_time = end_time - start_time;
	std::pair<bool, double> result(solution, running_time.count());
	return result;
}