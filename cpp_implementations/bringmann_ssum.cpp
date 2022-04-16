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
	
	int num_part = pow(2, ceil(log(layer / log_layer_delta) / log(2)));
	std::cout << log_layer_delta << " " << num_part << std::endl;
	
	std::vector<std::vector<int>> set_partitioned;

	std::vector<int> res;
	return res;
}


