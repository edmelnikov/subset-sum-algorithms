#include "fftw3.h"
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

#include "minkowski_sum.h"
#include "n_choose_k.h"



std::vector<std::vector<int>> rand_part(const std::vector<int>& set, int num_parts) {

	std::vector<std::vector<int>> partition(num_parts, std::vector<int>()); // partition of a set
																			
	// after both of these capacity changes back to zero?
	//for (auto subset : partition) {
	//	subset.reserve(set.size()); //  reserve enough space to avoid resizing during push_back()'s (maybe, it's a bit too much space?)
	//	
	//}
	//for (int i = 0; i < set.size(); i++) {
	//	partition.at(i).reserve(set.size());
	//}
	
	std::for_each(partition.begin(), partition.end(), [set](std::vector<int>& subset) { subset.reserve(set.size()); });

	int color;

	for (auto element : set) {
		color = rand() % num_parts; // assign a random color
		partition[color].push_back(element); // and push to the subset of that color
	}

	return partition;
}

std::vector<int> ColorCoding(const std::vector<int>& set, int target, int k, double delta) {
	
	int num_iter = ceil((log((1 / delta)) / log((double)4/3)));

	std::vector<std::vector<int>> set_partitioned;
	std::vector<int> mink_prod_united, mink_prod_united_old;
	// reserve space for all possible sums in mink_prod_united to avoid memory allocation
	// for that we will calculate number of combination 
	int num_sums = 0;
	for (int i = 0; i <= set.size(); i++) {  // set.size() could be replaced with target later 
		num_sums += n_choose_k(set.size(), i);
	}

	mink_prod_united.reserve(num_sums);
	mink_prod_united_old.reserve(num_sums);
	// std::cout << mink_prod_united.capacity() << std::endl;

	for (int i = 0; i < num_iter; i++) {
		set_partitioned = rand_part(set, k*k); // randomly partition into k^2 sets
		
		// and calculate minkowski sum for all the sets
		std::vector<int> minkowski_product = set_partitioned[0];
		for (int j = 1; j < k * k; j++) {  
			if (set_partitioned[j].size() > 0) {
				minkowski_product = minkowski_add(minkowski_product, set_partitioned[j], target); // NOTE: minkowsky sum always returns a sorted set
			}			
		}

		// unite minkowski sets
		mink_prod_united_old = mink_prod_united; // copy elements to the buffer "mink_prod_united_old" in order to record the result to "mink_prod_united"
		mink_prod_united.clear(); // clear to insert back the result of unity
		std::set_union(
			mink_prod_united_old.begin(),
			mink_prod_united_old.end(),
			minkowski_product.begin(),
			minkowski_product.end(),
			std::back_inserter(mink_prod_united)
		); 
		//for (auto el : minkowski_product) {
		//	std::cout << el << " ";
		//}
		//std::cout << std::endl;
	}
	
	return mink_prod_united;

}

std::vector<int> ColorCodingLayer(const std::vector<int>& set, int target, int layer, double delta) {
	double log_layer_delta = log((double)layer / delta) / log(2);
	if (layer < log_layer_delta) {
		return ColorCoding(set, target, layer, delta);
	}
	
	int num_part = pow(2, ceil(log(layer / log_layer_delta) / log(2))); // number of partitions
	// std::cout << log_layer_delta << " " << num_part << std::endl;
	std::vector<std::vector<int>> set_partitioned = rand_part(set, num_part);
	std::vector<std::vector<int>> set_partitioned_ssums(num_part, std::vector<int>()); // subset sums of a partitioned set

	// std::vector<int> subset_sums;

	for (int i = 0; i < num_part; i++) {
		//std::cout << "Colorcoding parameters: t=" << ceil(2 * 6 * log_layer_delta * target / layer);
		//std::cout << " k=" << ceil(6 * log_layer_delta) << " delta/layer=" << delta / layer << std::endl;

		set_partitioned_ssums[i] = ColorCoding(
			set_partitioned[i],
			ceil(2*6*log_layer_delta*target/layer),
			ceil(6*log_layer_delta),
			delta/layer
			);
	}

	// apply minkowski sum to all set_partitioned_ssums in a binary-tree-like way 
	std::vector<std::vector<int>> curr_round = set_partitioned_ssums; // copying
	
	/* calculate the highest number of ssums that can be sampled from the set */
	int num_sums = 0; 
	for (int i = 0; i <= set.size(); i++) {  
		num_sums += n_choose_k(set.size(), i);
	}	
	// std::vector<std::vector<int>> next_round(num_part, std::vector<int>());
	// std::for_each(next_round.begin(), next_round.end(), [num_sums](std::vector<int>& subset) { subset.reserve(num_sums); });
	std::vector<std::vector<int>> next_round;
	next_round.reserve(num_part);  // reserve to fit all sets without resizing

	for (int round = 0; round < log(num_part) / log(2); round++) { // iterate over levels of a tree from bottom to top (round = depth - level)
		for (int j = 0; j < curr_round.size(); j += 2) {
			if (j + 1 < curr_round.size()) {
				next_round.push_back(minkowski_add(curr_round[j], curr_round[j+1], pow(2, round + 1)*2*6*log_layer_delta*target/layer));
			}
			else {
				next_round.push_back(curr_round[j]);
			}
		}
		curr_round.clear();
		curr_round = next_round;
		next_round.clear();
	}

	//for (auto set : curr_round) {
	//	std::cout << "[";
	//	for (auto sum : set) {
	//		std::cout << sum << " ";
	//	}
	//	std::cout << "]" << std::endl;
	//}
	//std::cout << std::endl;

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

std::pair<bool, double> bringmann_ssum(const std::vector<int>& set, int target, double delta) {
	auto start_time = std::chrono::high_resolution_clock::now(); // time measurement

	/* Split into ⌈log n⌉ layers*/
	int num_layers = ceil(log(set.size()) / log(2)); // calculate number of layers

	// std::vector<std::vector<int>> set_layers(num_layers, std::vector<int>()); // create a vector with layers
	// std::for_each(set_layers.begin(), set_layers.end(), [set](std::vector<int>& layer) { layer.reserve(set.size()); });
	
	std::vector<int> layer_set;
	layer_set.reserve(set.size());
	int set_ind = 0;
	double right_bound = (double)target / pow(2, num_layers - 1); // right bound of current layer
	std::vector<int> layer_ssums;
	std::vector<int> all_ssums = {};

	for (int layer_ind = num_layers - 1; layer_ind >= 0; layer_ind--) { // iterate over layers backwards

		while (set[set_ind] <= right_bound && set_ind < set.size()) { // check if current set's item fits into the layer
			layer_set.push_back(set[set_ind]); // if it does, add it to current layer
			set_ind++; // and increment index of a set
		}
		right_bound = (double)target / pow(2, layer_ind - 1); // recalculate the right bound once all the items have been added to current layer

		//std::cout << "Color coding layer call with target " << target << ", layer " << pow(2, layer_ind + 1);
		//std::cout << ", error probability " << delta / num_layers << std::endl;
		//std::cout << "set " << std::endl;
		//for (auto el : layer_set) {
		//	std::cout << el << " ";
		//}
		//std::cout << std::endl;

		layer_ssums = ColorCodingLayer(layer_set, target, pow(2, layer_ind + 1), delta/num_layers); // calculate all subset sums of a layer
			
		all_ssums = minkowski_add(all_ssums, layer_ssums, target); // NOTE: minkowsky sum always returns a sorted set
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