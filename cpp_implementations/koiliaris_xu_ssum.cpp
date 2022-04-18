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

/* Function used in lemma 2.11 */
/*
	Given an element of a set B ⊆ [left_bound, ..., right_bound]×[0, ..., alpha], we compute the following function:
	f((i, j)) = (i - xj, j) 
*/
std::pair<int, int> function(const std::pair<int, int>& pair, int left_bound) {
	std::pair<int, int> res = {pair.first - pair.second*left_bound, pair.second};
	return res;
}

/* Inverse of the function above used in lemma 2.11 */
/*
	f^(-1)((i, j)) = (i + xj, j)
*/
std::pair<int, int> function_inverse(const std::pair<int, int>& pair, int left_bound) {
	std::pair<int, int> res = { pair.first + pair.second * left_bound, pair.second };
	return res;
}

/* Lemma 2.11 */
/*
	Given two (ssum, card)-sets (ssum_card_set1 and ssum_card_set2) with cardinality bound alpha based on disjoint sets 
	B, C ⊆ [left_bound, ..., right_bound] respectively, the following function computes (ssum, card)-set with cardinality 
	bound alpha for a union of sets B, C
*/
std::vector<std::pair<int, int>> ssum_card_set_union(const std::vector<std::pair<int, int>>& ssum_card_set1,
	const std::vector<std::pair<int, int>>& ssum_card_set2, int alpha, int left_bound, int right_bound) {

	/* Apply the function to the sets*/
	std::vector<std::pair<int, int>> ssum_card_set1_f; // here we will store the function outputs for all elements of a set
	ssum_card_set1_f.reserve(ssum_card_set1.size());
	std::vector<std::pair<int, int>> ssum_card_set2_f;
	ssum_card_set2_f.reserve(ssum_card_set2.size());

	for (auto el : ssum_card_set1) {
		ssum_card_set1_f.push_back(function(el, left_bound));
	}
	for (auto el : ssum_card_set2) {
		ssum_card_set2_f.push_back(function(el, left_bound));
	}

	/* Apply minkowski sum to the sets of function outputs */
	std::vector<std::pair<int, int>> min_add = minkowski_add_2d(ssum_card_set1_f, ssum_card_set2_f);

	/* exclude out-of-bounds values ((X ⊕ Y ) ∩ ([0 : ℓα] × [0 : α]))*/
	min_add.erase(
		std::remove_if(
			min_add.begin(),
			min_add.end(),
			[alpha, left_bound, right_bound](const std::pair<int, int>& el) { return el.first > alpha*(right_bound - left_bound) || el.second > alpha; }),
		min_add.end()
	);

	/* Apply the inverse function to the minkowski sum bounded product*/
	std::vector<std::pair<int, int>> min_add_f_inv;
	min_add_f_inv.reserve(min_add.size());
	for (auto el : min_add) {
		min_add_f_inv.push_back(function_inverse(el, left_bound));
	}

	return min_add_f_inv;
}

/* Lemma 2.12 */
/*
	Given a set and bound alpha the following function computes a (ssum, card)-set with cardinality bound alpha for the set
*/
std::vector<std::pair<int, int>> compute_ssum_card_set(const std::vector<int>& set, int alpha) {
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
	std::vector<std::pair<int, int>> left_ssum_card_set = compute_ssum_card_set(left_subset, alpha);
	std::vector<std::pair<int, int>> right_ssum_card_set = compute_ssum_card_set(right_subset, alpha);
	
	/* Once calculated left and right (ssum, card)-sets, combine them using ssum_card_set_union() (lemma 2.11) */
	return ssum_card_set_union(left_ssum_card_set, right_ssum_card_set, alpha, left_bound, right_bound);
}

/* Lemma 2.13 */
/*
	Given a set set S ⊆ [left_bound, ..., right_bound] (i. e., a set IN RANGE from left_bound to right_bound, hence the name of the function),
	the following function computes all subset sums less than or equal to ssum_bound using (ssum, card)-set
	TODO: get rid of sorting either by utilizing std::set's or by changing compute_ssum_card_set() function
*/
std::vector<int> compute_ssums_in_range(const std::vector<int>& set, int ssum_bound) {

	int alpha = floor((double)ssum_bound/set.front()); // calculate cardinality bound 
	std::vector<std::pair<int, int>> ssum_card_set = compute_ssum_card_set(set, alpha); // calculate a (ssum, card)-set
	
	std::vector<int> ssums;
	ssums.reserve(ssum_card_set.size());
	for (auto el : ssum_card_set) {
		if (el.first <= ssum_bound) {
			ssums.push_back(el.first);
		}
	}
	std::sort(ssums.begin(), ssums.end());

	return ssums;
}

/* Theorem 2.7 */
/*
	Given a set with total sum σ, the following function computes all of its subset sums in O(σ log(σ)log(n))
*/
std::vector<int> compute_ssums(const std::vector<int>& set, int ssum_bound) {
	if (set.size() == 1) return std::vector<int>({set.front()}); // if a set contains only one element

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

}

// std::vector<int> compute_ssums()
