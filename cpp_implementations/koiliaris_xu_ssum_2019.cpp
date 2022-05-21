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



std::vector<std::pair<int, int>> compute_ssum_card_set(const std::vector<int>& set, int bound=-1) {
	if (set.size() == 0) return std::vector<std::pair<int, int>>({ {0, 0} });
	if (set.size() == 1) return std::vector<std::pair<int, int>>({{set[0], 1}}); // if a set contains only one element
	int median_ind = set.size() / 2; // index of a set median
	/* Bounds that a set belongs to */
	int left_bound = set.front();
	int right_bound = set.back();
	
	/* Construct left and right subsets*/
	std::vector<int>::const_iterator subset_start;
	std::vector<int>::const_iterator subset_end;
		/* left */
	subset_start = set.begin();
	subset_end = set.begin() + median_ind;
	std::vector<int> left_subset(subset_start, subset_end);
		/* right */
	subset_start = set.begin() + median_ind;
	subset_end = set.end();
	std::vector<int> right_subset(subset_start, subset_end);

	/* Calculate (ssum, card)-sets for left and right subsets */
	std::vector<std::pair<int, int>> left_ssum_card_set = compute_ssum_card_set(left_subset, bound);
	std::vector<std::pair<int, int>> right_ssum_card_set = compute_ssum_card_set(right_subset, bound);
	
	// std::cout << "Lengths " << left_ssum_card_set.size() << " " << right_ssum_card_set.size() << std::endl;
	return minkowski_sum_fft2d(left_ssum_card_set, right_ssum_card_set, bound);
}

// lemma 3.3
std::vector<int> compute_ssums_congr_class(const std::vector<int>& set, int divisor, int target) {
	if (set.size() == 0 || set.size() == 1) return set;
	std::vector<int> whole_quots; // whole quotients
	whole_quots.reserve(set.size());

	for (auto el : set) {
		whole_quots.push_back(el / divisor);
	}

	std::vector<std::pair<int, int>> ssum_card_set = compute_ssum_card_set(whole_quots, target / divisor);

	int remainder = set.front() % divisor;
	std::vector<int> ssums;
	for (auto el : ssum_card_set) {
		ssums.push_back(el.first*divisor + el.second*remainder);
	}
	return ssums;
}


// theorem 3.4
std::pair<bool, double> koiliaris_xu_ssum(const std::vector<int>& set, int target) {
	auto start_time = std::chrono::high_resolution_clock::now(); // time measurement
	double set_size = set.size();
	// int num_part = floor(sqrt(set_size * (log(set_size)/ log(2))));
	int num_part = floor(sqrt(set_size * log(set_size)));

	/* Partition the set */
	std::vector<std::vector<int>> set_part(num_part, std::vector<int>()); // set partitione
	int remainder;
	for (auto el : set) {
		remainder = el % num_part;
		set_part[remainder].push_back(el);
	}

	std::vector<int> ssums = {};
	for (int congr_class = 0; congr_class < num_part; congr_class++) {
		ssums = minkowski_sum_fft(
			ssums,
			compute_ssums_congr_class(set_part[congr_class], num_part, target),
			target
			);
	}
	/*
	for (auto el : ssums) {
		std::cout << el << ", ";
	}
	std::cout << std::endl;*/



	auto end_time = std::chrono::high_resolution_clock::now(); // time measurement
	std::chrono::duration<double, std::milli> running_time = end_time - start_time;
	std::pair<bool, double> result(
		(target == ssums.back()),
		running_time.count());

	return result;
}