#pragma once
# include <vector>

std::pair<int, int> function(const std::pair<int, int>& pair, int left_bound);
std::pair<int, int> function_inverse(const std::pair<int, int>& pair, int left_bound);
std::vector<std::pair<int, int>> ssum_card_set_union(const std::vector<std::pair<int, int>>& ssum_card_set1,
	const std::vector<std::pair<int, int>>& ssum_card_set2, int alpha, int left_bound, int right_bound);
//std::vector<std::pair<int, int>> compute_ssum_card_set(const std::vector<int>& set, int alpha);
std::vector<std::pair<int, int>> compute_ssum_card_set(const std::vector<int>& set, int bound=-1);
std::vector<int> compute_ssums_in_range(const std::vector<int>& set, int ssum_bound);
std::vector<int> compute_ssums(const std::vector<int>& set, int ssum_bound);
std::vector<std::vector<int>> partition(const std::vector<int>& set, double r0);
std::pair<bool, double> koiliaris_xu_ssum(const std::vector<int>& set, int target, double r0 = -1);